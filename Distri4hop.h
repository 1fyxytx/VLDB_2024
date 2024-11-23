#include "../utils/communication.h"
#include "../blogel/BVertex.h"
#include "../blogel/Block.h"
#include "../blogel/BWorker.h"
#include "../blogel/BGlobal.h"
#include "../blogel/BType.h"
#include <queue>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <ext/hash_map>
#include <string.h>
#include <map>
#include <numeric>
#include <utility>
#include <time.h>
#include <deque>
#include <set>
#include <algorithm>
#include <bitset>
#include <vector>
#include <omp.h>

// 核心点：并行化树分解，如何确认顶点之间的互不影响性
// 结合分区和树分解的性质，进一步加速分解过程
#define MAXINT ((unsigned) 4294967295)

using namespace std;

struct EnumValue {
	vector<int> Local, Edge; 

    void set_Edge(int adj_id, int id1){

		if (v2part[id1] != v2part[adj_id])
			Edge.push_back(adj_id);
		else
			Local.push_back(adj_id);
    }

	void empty(){
		vector<int>().swap(Local);
		vector<int>().swap(Edge);
	}

    friend ibinstream & operator<<(ibinstream & m, const EnumValue & v){
		m<<v.Local;
		m<<v.Edge;

    	return m;
    }


    friend obinstream & operator>>(obinstream & m, EnumValue & v){
		m>>v.Local;
    	m>>v.Edge;
    	
		return m;
    }
};


struct MsgInf {
	vector<treeinf> Neig;

	MsgInf(){}

    friend ibinstream& operator<<(ibinstream& m, const MsgInf& idm)
    {
        return m;
    }

    friend obinstream& operator>>(obinstream& m, MsgInf& idm)
    {
        return m;
    }
};


class LCRVertex : public BVertex<VertexID, EnumValue, MsgInf>{
public:
	virtual void compute(MessageContainer& messages){}

	// === 构建 dist ==

	


	// =======边界索引构建策略=========、
	void Bound_Compute(MessageContainer& messages){ // 按道理讲，这个过程不需要接收其他分区顶点的信息

		vote_to_halt();
	}

};


class LCRBlock : public Block<char, LCRVertex, MsgInf> {
public:

	void Boundary_compute(VertexContainer &vertexes){ // 边界图计算,不接收信息

		vote_to_halt();
	}



	virtual void compute(MessageContainer &messages, VertexContainer &vertexes){}

};


class LCRBlockWorker : public BWorker<LCRBlock>{
public:
	unsigned MAXDIS, MAXMOV, MASK;
	vector<unsigned> *label, *flags;
	vector<vector<unsigned> > con, conB;
	bool *usd_bp;
	bool* is_indep;
	BPLabel *label_bp;
    int BoundNum = 0;

	virtual LCRVertex *toVertex(char *line){
		LCRVertex *v;
		return v;
	}


	void construct_bp_label(){
        label_bp = new BPLabel[totalV];
        usd_bp = new bool[totalV];
        memset( usd_bp, 0, sizeof(bool) * totalV );
        vector<int> v_vs[N_ROOTS];

        // ===============================  //
        int r = 0;
        for (int i_bpspt = 0; i_bpspt < N_ROOTS; ++i_bpspt) {
            while (r < totalV && usd_bp[r]) ++r;
            if (r == totalV) {
                for (int v = 0; v < totalV; ++v) 
                    label_bp[v].bpspt_d[i_bpspt] = MAXD;
                continue;
            }
            usd_bp[r] = true;
            v_vs[i_bpspt].push_back(r);
            int ns = 0;
            for (int i = 0; i < con[r].size(); ++i) {
                int v = con[r][i];
                if (!usd_bp[v]) {
                    usd_bp[v] = true;
                    v_vs[i_bpspt].push_back(v);
                    if (++ns == 64) break;
                }
            }
        }

        // ===============================  //
        int n_threads = 1;
        #pragma omp parallel
        {
            if(omp_get_thread_num() == 0) n_threads = omp_get_num_threads();
        }
        if( n_threads > MAX_BP_THREADS ) omp_set_num_threads(MAX_BP_THREADS);

        #pragma omp parallel
        {
            int pid = omp_get_thread_num(), np = omp_get_num_threads();
            // if( pid == 0 ) 
            //     printf( "n_threads_bp = %d\n", np );
            vector<uint8_t> tmp_d(totalV);
            vector<pair<uint64_t, uint64_t> > tmp_s(totalV);
            vector<int> que(totalV);
            vector<pair<int, int> > child_es(totalEdges*2);

            for (int i_bpspt = pid; i_bpspt < N_ROOTS; i_bpspt+=np) {
                // printf( "[%d]", i_bpspt );

                if( v_vs[i_bpspt].size() == 0 ) continue;
                fill(tmp_d.begin(), tmp_d.end(), MAXD);
                fill(tmp_s.begin(), tmp_s.end(), make_pair(0, 0));

                r = v_vs[i_bpspt][0];
                int que_t0 = 0, que_t1 = 0, que_h = 0;
                que[que_h++] = r;
                tmp_d[r] = 0;
                que_t1 = que_h;

                for( size_t i = 1; i < v_vs[i_bpspt].size(); ++i) {
                    int v = v_vs[i_bpspt][i];
                    que[que_h++] = v;
                    tmp_d[v] = 1;
                    tmp_s[v].first = 1ULL << (i-1);
                }

                for (int d = 0; que_t0 < que_h; ++d) {
                    int num_child_es = 0;

                    for (int que_i = que_t0; que_i < que_t1; ++que_i) {
                        int v = que[que_i];

                        for (int i = 0; i < con[v].size(); ++i) {
                            int tv = con[v][i];
                            int td = d + 1;

                            if (d == tmp_d[tv]) {
                                if (v < tv) {
                                    tmp_s[v].second |= tmp_s[tv].first;
                                    tmp_s[tv].second |= tmp_s[v].first;
                                }
                            } else if( d < tmp_d[tv]) {
                                if (tmp_d[tv] == MAXD) {
                                    que[que_h++] = tv;
                                    tmp_d[tv] = td;
                                }
                                child_es[num_child_es].first  = v;
                                child_es[num_child_es].second = tv;
                                ++num_child_es;
                          }
                        }
                    }

                    for (int i = 0; i < num_child_es; ++i) {
                        int v = child_es[i].first, c = child_es[i].second;
                        tmp_s[c].first  |= tmp_s[v].first;
                        tmp_s[c].second |= tmp_s[v].second;
                    }

                    que_t0 = que_t1;
                    que_t1 = que_h;
                }

                for (int v = 0; v < totalV; ++v) {
                    label_bp[v].bpspt_d[i_bpspt] = tmp_d[v];
                    label_bp[v].bpspt_s[i_bpspt][0] = tmp_s[v].first;
                    label_bp[v].bpspt_s[i_bpspt][1] = tmp_s[v].second & ~tmp_s[v].first;
                }
            }		
        }
    }


    bool prune_by_bp(int u, int v, int d) {
        BPLabel &idx_u = label_bp[u], &idx_v = label_bp[v];
        for (int i = 0; i < N_ROOTS; ++i) {
            int td = idx_u.bpspt_d[i] + idx_v.bpspt_d[i];
            if (td - 2 <= d)
                td += (idx_u.bpspt_s[i][0] & idx_v.bpspt_s[i][0]) ? -2 :
                        ((idx_u.bpspt_s[i][0] & idx_v.bpspt_s[i][1]) | (idx_u.bpspt_s[i][1] & idx_v.bpspt_s[i][0])) ? -1 : 0;
            if (td <= d) return true;
        }
        return false;
    }


	virtual void blockInit(VertexContainer &, BlockContainer &blocks){
		MAXDIS = 2; MAXMOV = 1;
		while( MAXINT / (totalV * 2) >= MAXDIS ) {
			MAXDIS *= 2;
			++MAXMOV;
		}
		MASK = MAXDIS - 1;

        v2p.resize(totalV, -1);
		p2v.resize(totalV, -1);

        for (int pos=0; pos<vertexes.size(); ++pos){
            LCRVertex * v = vertexes[pos];
            
			int vid = v->id, flg = v2degree[vid];
			
            if (v->value().Edge.size() > 0) 
				BoundNum += 1;

			sortList.push_back(make_pair(flg, -vid)); // pos就是id
		}

        sort(sortList.rbegin(), sortList.rend());
		for (int i=0; i<sortList.size(); ++i){
			v2p[-sortList[i].second] = i;
			p2v[i] = -sortList[i].second; // new 2 old
		}

        con.resize(sortList.size());
		conB.resize(totalV);

        for (int i=0; i<sortList.size(); ++i){
            int ovid = p2v[i];
			// if (_my_rank == 1)
			// 	cout<<"old:  "<<ovid<<"   new:  "<<i<<endl;

            vector<int>& local = vertexes[vert2place[ovid]]->value().Local;
            for( int p = 0; p < local.size(); ++p ) {
                int jj = local[p]; 
                int j = v2p[jj];

                con[j].push_back(i);
		    }
			
			vector<int>& edge = vertexes[vert2place[ovid]]->value().Edge;
            for( int p = 0; p < edge.size(); ++p ) {
                int jj = edge[p]; 
				unsigned val = (jj << MAXMOV) | 1;

                conB[ovid].push_back(val); // 存的是原始id信息
		    }

			vertexes[vert2place[ovid]]->value().empty();
        }

		vector<pair<int, int>>().swap(sortList);
		vector<LCRVertex *>().swap(vertexes);
	}


	void AddVertex(char *line, int va){
		// 三个元素  V_A, V_B, Label
		int vb, pa, pb;
		// pa = v2part[va];
		LCRVertex* v = new LCRVertex;
		v->id = va; v->bid = 0;
		load_vertex(v);
		vert2place[va] = vertexes.size()-1;
		char* s1 = strtok(line," ");

		while(s1){
			totalEdges += 1;
			vb = atoi(s1)-1;
			// pb = v2part[vb];
			v->value().set_Edge(vb, va);			
			s1=strtok(NULL," ");
		}
		
        int flg = v->value().Edge.size()>0? 1:0;
        int vals = flg*10000000 + v->value().Edge.size() + 
                   v->value().Local.size();
		
		v2degree[va] = vals;
	}


	bool can_update(int v, int dis, char *nowdis) {
		for( int i = 0; i < (int) label[v].size(); ++i ) {
			int w = label[v][i]>>MAXMOV, d = label[v][i]&MASK;
			if( nowdis[w] >= 0 && nowdis[w] + d <= dis ) return false;
		}
		return true;
	}


	void Part2hop(){
		totalV = con.size();
		bool consider_indep = false;
        construct_bp_label();

		omp_set_num_threads(threads);
		double t = omp_get_wtime();

		int **pos = new int*[totalV];
		for( int i = 0; i < totalV; ++i ) 
			pos[i] = new int[MAXDIS];

		
		label = new vector<unsigned>[totalV]; // unsigned 整合了 id+dis
		flags = new vector<unsigned>[totalV]; 

		for(int i = 0; i < totalV; ++i){
			label[i].push_back((((unsigned)i)<<MAXMOV) | 0);
			flags[i].push_back(0); // 0表示这个index是必要的
		
			for( int j = 0; j < con[i].size() && con[i][j] < i; ++j ){
				label[i].push_back((((unsigned) con[i][j]) << MAXMOV) | 1);
				flags[i].push_back(0); // 0-1的标签都是直接相连的
			}

			pos[i][0] = 1, pos[i][1] = (int) label[i].size();
		}

		int dis = 2;
		for( long long cnt = 1; cnt && dis <= MAXDIS; ++dis ){
			
			cnt = 0;
			vector<unsigned> *label_new = new vector<unsigned>[totalV];
			vector<unsigned> *flg_new = new vector<unsigned>[totalV];
			#pragma omp parallel
			{
				int pid = omp_get_thread_num(), np = omp_get_num_threads();
				long long local_cnt = 0;
				unsigned char *used = new unsigned char[totalV/8+1];
				memset( used, 0, sizeof(unsigned char) * (totalV/8+1) );
				vector<int> cand, candFg, flgs(totalV);
				
				char *nowdis = new char[totalV];
				memset( nowdis, -1, sizeof(char) * totalV);
	
				for( int u = pid; u < totalV; u += np ){
					if (usd_bp[u]) continue;
					cand.clear();
					candFg.clear();
					
					for( int i = 0; i < con[u].size(); ++i ){
						int w = con[u][i];

						for( int j = pos[w][dis-2]; j < pos[w][dis-1]; ++j ) {
							int v = label[w][j] >> MAXMOV;
							if( v >= u ) break;

							int ff = flags[w][j];

							if (w < BoundNum) ff = 1; // 由另外的边界点传输而来
							
							if (ff == 1 and flgs[v] == 0){
								flgs[v] = 1;
								candFg.push_back(v);
							} 

							if( !(used[v/8]&(1<<(v%8))) ) {
								used[v/8] |= (1<<(v%8)), cand.push_back(v);
							}
					
						}
					}

					int n_cand = 0;
					for( int i = 0; i < (int) label[u].size(); ++i ) 
						nowdis[label[u][i] >> MAXMOV] = label[u][i] & MASK;
					
					for( int i = 0; i < (int) cand.size(); ++i ) {
						used[cand[i]/8] = 0;
						if( !prune_by_bp(u, cand[i], dis) )
							if( can_update(cand[i], dis, nowdis) ) 
								cand[n_cand++] = cand[i]; 
					}

					cand.resize(n_cand);
					sort(cand.begin(), cand.end());
					for( int i = 0; i < (int) cand.size(); ++i ) {
						label_new[u].push_back((((unsigned)cand[i])<<MAXMOV) | (unsigned) dis), ++local_cnt;
						flg_new[u].push_back(flgs[cand[i]]);
					}
					for( int i = 0; i < (int) label[u].size(); ++i ) 
						nowdis[label[u][i]>>MAXMOV] = -1;

					for (int vid: candFg) flgs[vid] = 0;
				}
				#pragma omp critical
				{
					cnt += local_cnt;
				}
				delete[] used; delete[] nowdis;
			}

			#	pragma omp parallel
			{
				int pid = omp_get_thread_num(), np = omp_get_num_threads();
				for( int u = pid; u < totalV; u += np ){
					label[u].insert(label[u].end(), label_new[u].begin(), label_new[u].end());
					flags[u].insert(flags[u].end(), flg_new[u].begin(), flg_new[u].end());

					vector<unsigned>(label[u]).swap(label[u]);
					vector<unsigned>(flags[u]).swap(flags[u]);

					vector<unsigned>().swap(label_new[u]);
					vector<unsigned>().swap(flg_new[u]);

					pos[u][dis] = (int) label[u].size();
				}
			}
			
			if (_my_rank == 0)
			cout<<"Distance: "<<dis<<"   Cnt: "<<cnt<<endl;
			delete[] label_new,flg_new;
		}

		delete[] pos;

		for(int i=BoundNum; i<totalV; ++i){
			Wavg += label[i].size();
			W1 += label[i].size();
		}

		

	}


	void Core2hop(){
		float ct = 0.02;
		totalV = con.size();
		int Tcnt = totalV*ct;

		bool consider_indep = false;

		omp_set_num_threads(threads);
		double t = omp_get_wtime();
		vector<int>* pos = new vector<int>[totalV];

		label = new vector<unsigned>[totalV]; // unsigned 整合了 id+dis
		
		for(int i = 0; i < totalV; ++i){
			if (i%_num_workers == _my_rank or i < Tcnt){
				label[i].push_back((((unsigned)i)<<MAXMOV) | 0);
				pos[i].push_back(1);
			}else{
				pos[i].push_back(0);
			}
		}

		int dis = 1;
		for( long long cnt = 1; cnt && dis <= MAXDIS; ++dis ){
			long long cnt2 = 0;
			cnt = 0;
			vector<unsigned> *label_new = new vector<unsigned>[totalV];
			#pragma omp parallel
			{
				int pid = omp_get_thread_num(), np = omp_get_num_threads();
				long long local_cnt = 0, local_candidate = 0;
				unsigned char *used = new unsigned char[totalV/8+1];
				memset( used, 0, sizeof(unsigned char) * (totalV/8+1) );
				vector<int> cand;
				char *nowdis = new char[totalV];
				memset( nowdis, -1, sizeof(char) * totalV);
	
				for( int u = pid; u < totalV; u += np ){
					// if (usd_bp[u]) continue;
					cand.clear();

					for( int i = 0; i < con[u].size(); ++i ){
						int w = con[u][i] >> MAXMOV, 
						    d = con[u][i] &  MASK;

						// cout<<u<<"  "<<w<<" "<<d<<endl;
						if (d == 0) continue;

						if (d > dis){
							local_candidate += 1;
							continue;
						}
						
						int lastVal = pos[w][pos[w].size()-1];

						if (lastVal - pos[w][dis-d] > 0){ // == 这些是还没有被判定的标签 ==
							local_candidate += 1;
						}

						for( int j = d==dis? 0:pos[w][dis-d-1]; j < pos[w][dis-d]; ++j ) {
							// if (j >= label[w].size())
							// cout<<"d: "<<d<<"  s1: "<<j<<"  s2: "<<pos[w][dis-d]<<"  label: "<<label[w].size()<<endl;
							int v = label[w][j] >> MAXMOV;
							if( v >= u ) break;

							if( !(used[v/8]&(1<<(v%8))) ) 
								used[v/8] |= (1<<(v%8)), cand.push_back(v);
						}
					}

					int n_cand = 0;
					for( int i = 0; i < (int) label[u].size(); ++i ) 
						nowdis[label[u][i] >> MAXMOV] = label[u][i] & MASK;
					
					for( int i = 0; i < (int) cand.size(); ++i ) {
						used[cand[i]/8] = 0;
						// if( !prune_by_bp(u, cand[i], dis) )
						if( can_update(cand[i], dis, nowdis) ) 
							cand[n_cand++] = cand[i]; 
					}

					cand.resize(n_cand);
					sort(cand.begin(), cand.end());
					for( int i = 0; i < (int) cand.size(); ++i ) {
						label_new[u].push_back((((unsigned)cand[i])<<MAXMOV) | (unsigned) dis), ++local_cnt;
					}
					for( int i = 0; i < (int) label[u].size(); ++i ) 
						nowdis[label[u][i]>>MAXMOV] = -1;
				}
				#pragma omp critical
				{
					cnt += local_cnt;
					cnt2 += local_candidate;
				}
				delete[] used; delete[] nowdis;
			}

			#	pragma omp parallel
			{
				int pid = omp_get_thread_num(), np = omp_get_num_threads();
				for( int u = pid; u < totalV; u += np ){
					label[u].insert(label[u].end(), label_new[u].begin(), label_new[u].end());
					vector<unsigned>(label[u]).swap(label[u]);
					vector<unsigned>().swap(label_new[u]);
					pos[u].push_back((int) label[u].size());
				}
			}

			if (_my_rank == 0)
				cout<<"Distance: "<<dis<<"   Cnt: "<<cnt<<"  "<<cnt2<<endl;
			
			cnt += cnt2;
			delete[] label_new;
		}

		// vector<vector<unsigned> >().swap(con);
		// for(int i = 0; i < totalV; ++i){
		// 	label[i].push_back((((unsigned)i)<<MAXMOV) | 0);
		// }
	}


	void LabelPrune(){ // == 有些标签删除之后，对于距离是不会有影响的，可并行 ==
		long long total=0, cnt=0;
		for (int i=0; i<BoundNum; ++i){
			vector<unsigned>& lbs = flags[i];
			// ==== revert search ===
			for (int i=0; i<lbs.size(); ++i){
				if (lbs[i] == 0) cnt += 1;
			}
			total += lbs.size();
		}

		long long T1 = all_sum_LL(total), T2 = all_sum_LL(cnt);
		if (_my_rank == 0) cout<<"necessary:  "<<float(T2)/T1<<endl;
	}


	void BoundGraph(){
		for (int i=0; i<BoundNum; ++i){
			int oid = p2v[i];
			
			vector<unsigned>& lbs = label[i];
			vector<unsigned>& fls = flags[i];

			for (int j=0; j<lbs.size(); ++j){
				
				if (fls[j] == 1) continue; // 不影响距离

				unsigned vid = lbs[j] >> MAXMOV,
						 dis = lbs[j] & MASK;
				int ovid = p2v[vid];
				
				if (oid == ovid or dis == 0) continue;

				// dis = 1;
				unsigned elem1 = ovid << MAXMOV | dis,
						 elem2 = oid << MAXMOV | dis;
				// pair<int, int> elem1(ovid, dis), elem2(oid, dis);

				conB[oid].push_back(elem1);
				conB[ovid].push_back(elem2);
			}

			vector<unsigned>().swap(lbs);
			vector<unsigned>().swap(fls);
		}

		delete[] label;
		delete[] flags;

		long long edges = 0, TotalE = 0;
		for (int i=0; i<conB.size(); ++i)
			edges += conB[i].size();
		
		TotalE = all_sum_LL(edges);
		
		if (_my_rank == 0){
			cout<<"Bound V: "<<BdV<<"  E: "<<(long long)(TotalE/2)<<endl;
		}


		vector<vector<vector<unsigned > > > CB(_num_workers);
		

		// == 收集 CopyNum 顶点的标签信息  ==
		int start = 0, batch = 10000;
		if (start+batch > BoundNum)
			batch = BoundNum-start;

		while (1){
			for (int i=start; i<start+batch; ++i){
				int vid = p2v[i];
				conB[vid].push_back(vid);
				for (int j=0; j<_num_workers; ++j)
					CB[j].push_back(conB[vid]);
			}
			
			all_to_all(CB);
			worker_barrier();

			for(int i=0; i<_num_workers; ++i){
				for (int j=0; j<CB[i].size(); ++j){
					vector<unsigned>& inf = CB[i][j];
					int vid = inf[inf.size()-1];
					inf.pop_back();
					swap(conB[vid], inf);
					vector<unsigned>(conB[vid]).swap(conB[vid]);
				}
			}
			
			for (int i=0; i<_num_workers; ++i)
				vector<vector<unsigned > >().swap(CB[i]);

			vector<vector<vector<unsigned > > >(CB).swap(CB);

			start += batch;

			int cntt = all_sum(start);
			if (_my_rank == 0)
				cout<<"Cnt: "<<cntt<<"  "<<BdV<<endl;

			if (cntt == BdV) break;

			if (start+batch > BoundNum)
				batch = BoundNum-start;
		}

		totalV = ALLVertices;

		// == conB 格式转化 ==
		v2p.resize(totalV, -1);
		p2v.resize(totalV, -1);
		vector<pair<int, int>>().swap(sortList);

		for(int i=0; i<conB.size(); ++i){
			if (conB[i].size() == 0) 
				continue;
			sortList.push_back(make_pair(v2degree[i], -i));
		}

		sort(sortList.rbegin(), sortList.rend());
		for (int i=0; i<sortList.size(); ++i){
			v2p[-sortList[i].second] = i;
			p2v[i] = -sortList[i].second;
		}
		vector<vector<unsigned>>().swap(con);
		con.resize(sortList.size());

		for (int i=0; i<sortList.size(); ++i){
			int ovid = p2v[i];
			// cout<<"old:  "<<ovid<<"   new:  "<<i<<endl;
			vector<unsigned>& local = conB[ovid];
			for( int p = 0; p < local.size(); ++p ) {
				unsigned jj = local[p] >> MAXMOV,
						 dd = local[p] & MASK;

				unsigned j = v2p[jj], val = ((unsigned)i<<MAXMOV) | dd;
				
				// if (j >= con.size()) cout<<"error!"<<endl;
				con[j].push_back(val);
			}
		}
		vector<vector<unsigned>>().swap(conB);
		vector<pair<int, int>>().swap(sortList);
	
		for (int i=0; i<con.size(); ++i){
			vector<unsigned>(con[i]).swap(con[i]);
		}

		// // string new_filename = "/home/hnu/Disk0/zyy_dataset/uk2002/graph.txt";
		// // const char *file2 = new_filename.c_str();
		// // fstream outfileX;
		// // outfileX.open(file2, ios::out);
		// // for (int i=0; i<con.size(); ++i){
		// // 	for (int j=0; j<con[i].size(); ++j){
		// // 		int id = con[i][j] >> MAXMOV;
		// // 		outfileX<<to_String(i)<<" "<<to_String(id)<<endl;
		// // 	}
		// // }



		worker_barrier();
		// cout<<"End!! "<<con.size()<<endl;
	}


	void IndexPrint(){
		for (int i=0; i<con.size(); ++i){

			cout<<" id:  "<<p2v[i]<<endl;
			cout<<"$";
			for (int j=0; j<label[i].size(); ++j){
				unsigned val = label[i][j];
				unsigned vid = p2v[val>>MAXMOV], dis = val&MASK;
				cout<<"\\{v_"<<vid<<","<<dis<<"\\},";
			}
			cout<<"$"<<endl;
		}
	}


	void DHI_Store(string new_filename){
		FILE *fout = fopen( (new_filename).c_str(), "wb" );
		int cntt = con.size()-BoundNum;
		fwrite(&cntt, sizeof(int), 1, fout);

		for (int i=BoundNum; i<con.size(); ++i){
			int len = (int) label[i].size();
			fwrite(&len, sizeof(int), 1, fout);
		}

		for (int i=BoundNum; i<con.size(); ++i){
			int len = (int) label[i].size();
			unsigned *s = new unsigned[len]; // 第一个存放自身id

			for(int j = 0; j < len; ++j){
				unsigned vid = label[i][j]>>MAXMOV, dis = label[i][j]&MASK;
				unsigned vv = unsigned(p2v[vid]);
				// if (_my_rank == 0)
				// cout<<vid<<"  "<<vv<<"  "<<dis<<endl;
				unsigned val = vv << MAXMOV | dis;
				s[j] = val;
			}

			fwrite(s, sizeof(unsigned), len, fout);
			delete[] s;
		}

		fclose(fout);
	}


	void DHBound_Store(string new_filename){
		FILE *fout = fopen( (new_filename).c_str(), "wb" );
		int cntt = con.size();
		fwrite(&cntt, sizeof(int), 1, fout);

		for (int i=0; i<con.size(); ++i){
			int len = (int) label[i].size();
			fwrite(&len, sizeof(int), 1, fout);
		}

		for (int i=0; i<con.size(); ++i){
			int len = (int) label[i].size();
			unsigned *s = new unsigned[len]; // 第一个存放自身id

			for(int j = 0; j < len; ++j){
				unsigned vid = label[i][j]>>MAXMOV, dis = label[i][j]&MASK;
				unsigned vv = unsigned(p2v[vid]);
				// if (_my_rank == 0)
				// cout<<vid<<"  "<<vv<<"  "<<dis<<endl;
				unsigned val = vv << MAXMOV | dis;
				s[j] = val;
			}

			fwrite(s, sizeof(unsigned), len, fout);
			delete[] s;
		}

		fclose(fout);
	}


	void run_LCR(const WorkerParams& params){
		
		ifstream infile, infile1, infile2; // 先把分区文件读取
		
		threads = params.khop;

		const char *part_path = params.partition_path.c_str(); 
		infile.open(part_path);
		string s;
		if(!infile.is_open()){
			cout<<"No such file!"<<endl;
			exit(-1);
		}

		int nnn = 0;
		while(getline(infile, s)){
			char* part = new char[strlen(s.c_str())+1];
			strcpy(part, s.c_str());
			v2part.push_back( atoi(part) ); // 
			v2degree.push_back(0); 
			delete part;
			nnn += 1;
		}

		const char *filepath = params.input_path.c_str(); //"data//Amazon_New.txt";

		// ======================
		infile1.open(filepath);
		if(!infile1.is_open()){
			cout<<"No such file!"<<endl;
			exit(-1);
		}

		int xxx = 0;
		while(getline(infile1, s)){
			if (xxx > 0 and v2part[xxx-1] == _my_rank){
				char* strc = new char[strlen(s.c_str())+1];
				strcpy(strc, s.c_str());
				AddVertex(strc, xxx-1);
				delete strc;
			}
			xxx += 1;
		}

		vector<vector<int> > V2D(_num_workers);
		for(int i=0; i<_num_workers; ++i){
			if (i == _my_rank) 
				continue;
			V2D[i] = v2degree;
		}
		
		all_to_all(V2D);
		worker_barrier();

		for (int i=0; i<_num_workers; ++i){
			for (int j=0; j<V2D[i].size(); ++j){
				v2degree[j] += V2D[i][j]; 
			}
		}

		totalV = all_sum(vertexes.size());
		ALLVertices = all_sum(vertexes.size());
		long long totalE = all_sum_LL(totalEdges);

		blockInit(vertexes, blocks); 
		BdV = all_sum(BoundNum);

		if (_my_rank == 0){
			cout<<"Total V: "<<totalV<<"  E: "<<totalE/2<<endl;
			cout<<"Bound V: "<<BdV<<endl;
		}
		
		float t = omp_get_wtime();

		Part2hop();

		// string new_filename = params.output_path + "In_"+ to_String(_my_rank);
		// DHI_Store(new_filename);

		worker_barrier();

		BoundGraph();
		float ttt1 = omp_get_wtime()-t;

		Core2hop();

		// new_filename = params.output_path + "Bd_"+ to_String(_my_rank);
		// DHBound_Store(new_filename);
		
		float ttt = omp_get_wtime()-t;
		if (_my_rank == 0){
			cout<<"time: "<<ttt<<" s"<<endl;
			cout<<"Inner time: "<<ttt1<<" s"<<endl;

		}

		#	pragma omp parallel
		{
			long long wwag = 0;
			int pid = omp_get_thread_num(), np = omp_get_num_threads();

			for(int u = pid; u < totalV; u += np){
				for (int j=0; j<label[u].size(); ++j){
					int vid = label[u][j] >> MAXMOV;
					if (vid%_num_workers == _my_rank)
						wwag += 1;
				}
			}

			#pragma omp critical
			{
				Wavg += wwag;
			}
		}

		long long sss = all_sum_LL(Wavg);
		long long ss1 = all_sum_LL(W1);
		if (_my_rank == 0){
			cout<<"Memory size: "<<float(sss*4)/(1024*1024)<<" MB"<<endl;
			cout<<"Inner index size: "<<float(ss1*4)/(1024*1024)<<" MB"<<endl;
		}
	}
};



void SetLCR_Construction(string in_path, string partition_path, string out_path, int src, int dst, int khop){
	WorkerParams param;
	param.input_path=in_path;
	param.partition_path = partition_path;
	param.output_path=out_path;
	param.src = src;
	param.dst = dst;
	param.khop = khop;
	param.force_write=true;
	param.native_dispatcher=false;
	LCRBlockWorker worker;
	worker.set_compute_mode(LCRBlockWorker::VB_COMP);

	worker.run_LCR(param);
};














