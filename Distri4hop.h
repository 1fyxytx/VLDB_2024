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
#include <cmath>
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
	int tgt;
	unsigned inf; // vid + dis

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
	unsigned MAXDIS1, MAXMOV1, MASK1;

	vector<unsigned> *label, *flags, *root; // record the distance to landmarks
	vector<vector<unsigned> > con, conB;
	bool *usd_bp;
	BPLabel *label_bp;
    int BoundNum = 0, topk = 5;
	vector<pair<unsigned, unsigned> > CandV;
	int tstid = 0;
	long long LocaEdges = 0;
	vector<int> ActiveV;
	unsigned actV = 0;

	virtual LCRVertex *toVertex(char *line){
		LCRVertex *v;
		return v;
	}

	vector<int> Convert(unsigned temp)
	{
		vector<int> result;
		int i = 0;
		while (temp != 0)
		{
			if ((temp & 1) == 1)
			{
				result.push_back(i);
			}
			i++;
			temp = temp >> 1;
		}
		return result;
	}


	void merge(vector<pair<unsigned, unsigned>>& elem){
		sort(elem.begin(), elem.end());
		unsigned tgt = elem[0].first;
		int p = 0;

		for (int i=1; i<elem.size(); ++i){
			unsigned vid = elem[i].first;
			if (tgt == vid){
				elem[p].second |= elem[i].second;
			}else{
				p += 1;
				elem[p] = elem[i]; 
				tgt = vid;
			}
		}

		vector<pair<unsigned, unsigned> >(elem).swap(elem);
	}


	void construct_bp_label(){ // == vFlag 和 candidat ==
		int dis = 0;
		vector<vector<pair<unsigned, unsigned> > > Messages(_num_workers);
		vector<pair<unsigned, unsigned> > LocalInf;

		while (true){
			
			for (int i=0; i<CandV.size(); ++i){
				
				unsigned v1 = CandV[i].first, landV = CandV[i].second; // v2 是 topk个点的label的集合
				vector<int> TgtV = Convert(landV);

				int nowID = v2p[v1];
				
				vector<unsigned>& rlab = root[nowID];

				unsigned temp = 0;
				for (int j=0; j<TgtV.size(); ++j){
					int pla = TgtV[j];
					if (rlab[pla] == 100000){
						rlab[pla] = dis;
						unsigned val = pow(2,pla);
						temp |= val;						
					}
				}
				
				// == transfer to neighbors ==
				if (temp == 0) continue; // 不需要传标签

				for (int ii=0; ii<con[nowID].size(); ++ii){ // 同分区的顶点
					int vid = p2v[con[nowID][ii]]; // 
					LocalInf.push_back(make_pair(vid, temp));
				}

				int nowBndID = v2p_Global[v1];
				if (nowBndID == -1) continue;


				for (int ii=0; ii<conB[nowBndID].size(); ++ii){
					int vid = conB[nowBndID][ii] >> MAXMOV1,
						ovid = p2v_Global[vid];

					Messages[v2part[ovid]].push_back(make_pair(ovid, temp));
				}
			}

			// == 先完成一波整合，基于 | 运算 ==
			for (int i=0; i<_num_workers; ++i)
				if (Messages[i].size() > 1)
					merge(Messages[i]);

			all_to_all(Messages);
			worker_barrier();

			CandV.clear();	
			swap(LocalInf, CandV);

			for (int i=0; i<_num_workers; ++i){
				CandV.insert(CandV.end(), Messages[i].begin(), Messages[i].end());
				vector<pair<unsigned, unsigned> >().swap(Messages[i]);
			}

			if (CandV.size() > 1) 
				merge(CandV);
			
			long long cnt = all_sum_LL(CandV.size());
			worker_barrier();

			if (cnt == 0)
				break;

			if (_my_rank == 0)
				cout<<dis <<"  Cnt: "<<cnt<<endl;

			vector<vector<pair<unsigned, unsigned> > >().swap(Messages);
			Messages.resize(_num_workers);

			dis += 1;
		}
    }


    int prune_by_root(int u, int v, int d) {
		// 依托 root 的labels 进行剪枝
		int mind = 1000000;
		for (int i=0; i<topk; ++i){
			int dd = root[u][i]+root[v][i];
			if (mind > dd)
				mind = dd;
		}

		if (d >= mind) 
			return 0;
		else
			return 1;
    }


	virtual void blockInit(VertexContainer &, BlockContainer &blocks){
		// sortList, sortTotal confirm	
		for (int i=0; i<v2degree.size(); ++i){
			if (v2degree[i] > 100000000)
				sortTotal.push_back(make_pair(v2degree[i], -i)); // pos就是id

			if (v2part[i] == _my_rank)
				sortList.push_back(make_pair(v2degree[i], -i));
		}

		// cout<<_my_rank<<"  "<<sortTotal.size()<<endl;

		v2p.resize(totalV, -1);
		p2v.resize(totalV, -1);
		
		sort(sortList.rbegin(), sortList.rend());
		
		for (int i=0; i<sortList.size(); ++i){
			v2p[-sortList[i].second] = i;
			p2v[i] = -sortList[i].second; // new 2 old
		}

		v2p_Global.resize(totalV, -1); // id 不连续
		p2v_Global.resize(totalV, -1);
		
		sort(sortTotal.rbegin(), sortTotal.rend());
		
		for (int i=0; i<sortTotal.size(); ++i){
			v2p_Global[-sortTotal[i].second] = i;
			p2v_Global[i] = -sortTotal[i].second; // new 2 old
		}

		// ======================================
        con.resize(sortList.size());
		conB.resize(sortTotal.size());
		root = new vector<unsigned>[sortList.size()]; // 收集到landmark的distance值

		MAXDIS1 = 2; MAXMOV1 = 1;
		while( MAXINT / (BdV * 2) >= MAXDIS1 ) {
			MAXDIS1 *= 2;
			++MAXMOV1;
		}
		MASK1 = MAXDIS1 - 1;


		for (int i=0; i<sortList.size(); ++i){
            int ovid = p2v[i];
			root[i].resize(topk, 100000);

            vector<int>& local = vertexes[vert2place[ovid]]->value().Local;
            
			for( int p = 0; p < local.size(); ++p ) {
                int j = v2p[local[p]];

                con[j].push_back(i); // 内部存储的
		    }
			
			// ===================================
			vector<int>& edge = vertexes[vert2place[ovid]]->value().Edge;
			if (edge.size() > 0){
				int newPos = v2p_Global[ovid]; // 在边界图上的定位
				for( int p = 0; p < edge.size(); ++p ) {
					int newVid = v2p_Global[edge[p]]; // 在边界图上的定位

					unsigned elem = newVid << MAXMOV1 | 1;

					conB[newPos].push_back(elem); // 存的是原始id信息
				}
			}


			vertexes[vert2place[ovid]]->value().empty();
        }

		vector<pair<int, int>>().swap(sortList);
		vector<LCRVertex *>().swap(vertexes);


		// // ======================================
		for (int i=0; i<topk; ++i){
			int vid = - sortTotal[i].second;
			if (v2part[vid] == _my_rank){ 
				unsigned tgt = pow(2,i); // vid 当前激活点，i 目标点
				CandV.push_back(make_pair(vid, tgt));
			}
		}
		vector<pair<int, int>>().swap(sortTotal);
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
        int vals = flg*100000000 + v->value().Edge.size() + 
                   v->value().Local.size();
		
		if (flg == 1) BoundNum += 1;

		v2degree[va] = vals;
		LocaEdges += v->value().Local.size();
	}


	bool can_update(int v, int dis, char *nowdis) {
		for( int i = 0; i < (int) label[v].size(); ++i ) {
			int w = label[v][i]>>MAXMOV, d = label[v][i]&MASK;
			if( nowdis[w] >= 0 && nowdis[w] + d <= dis ) return false;
		}
		return true;
	}


	void Part2hop(){
		actV = *max_element(ActiveV.begin(), ActiveV.end());
		MAXDIS = 2; MAXMOV = 1;
		while( MAXINT / (actV * 2) >= MAXDIS ) {
			MAXDIS *= 2;
			++MAXMOV;
		}
		MASK = MAXDIS - 1;

		totalV = con.size();

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
				vector<int> cand;
				
				char *nowdis = new char[totalV];
				memset( nowdis, -1, sizeof(char) * totalV);
	
				for( int u = pid; u < totalV; u += np ){

					cand.clear();
					unordered_map<int, int> flg;					
					for( int i = 0; i < con[u].size(); ++i ){
						int w = con[u][i];

						for( int j = pos[w][dis-2]; j < pos[w][dis-1]; ++j ) {
							int v = label[w][j] >> MAXMOV;
							if( v >= u ) break;

							int ff = flags[w][j];

							if (w < BoundNum) ff = 1; // 由另外的边界点传输而来
							
							if (ff == 1) flg[v] = ff; // 由于领居顺序不一致，所以一旦有1标记，就不用改了

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

						int flg = prune_by_root(cand[i], u, dis);
						if (flg == 1)
							if( can_update(cand[i], dis, nowdis) ) 
								cand[n_cand++] = cand[i]; 
					}

					cand.resize(n_cand);
					sort(cand.begin(), cand.end());
					for( int i = 0; i < (int) cand.size(); ++i ) {
						label_new[u].push_back((((unsigned)cand[i])<<MAXMOV) | (unsigned) dis), ++local_cnt;
						flg_new[u].push_back(flg[cand[i]]);
					}
					for( int i = 0; i < (int) label[u].size(); ++i ) 
						nowdis[label[u][i]>>MAXMOV] = -1;
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
			
			if (_my_rank == tstid)
			cout<<"Distance: "<<dis<<"   Cnt: "<<cnt<<endl;

			delete[] label_new,flg_new;
		}

		delete[] pos;

		for(int i=BoundNum; i<totalV; ++i){
			Wavg += (label[i].size()+topk);
			W1 += (label[i].size()+topk);
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

		delete[] pos;
		// vector<vector<unsigned> >().swap(con);
		// for(int i = 0; i < totalV; ++i){
		// 	label[i].push_back((((unsigned)i)<<MAXMOV) | 0);
		// }
	}


	void BoundGraph(){

		for (int i=0; i<BoundNum; ++i){
			// obnd_id 是在边界图里面的 位置
			int oid = p2v[i], 
			    ob_id = v2p_Global[oid];
			
			vector<unsigned>& lbs = label[i];
			vector<unsigned>& fls = flags[i];

			for (int j=0; j<lbs.size(); ++j){
				
				if (fls[j] == 1) continue; // 不影响距离

				unsigned vid = lbs[j] >> MAXMOV,
						 dis = lbs[j] & MASK;

				int ovid = p2v[vid], ob_vid = v2p_Global[ovid];
				
				if (ob_id == ob_vid or dis == 0) continue;

				// ============= 
				unsigned elem1 = ob_vid << MAXMOV1 | dis,
					     elem2 =  ob_id << MAXMOV1 | dis;

				
				conB[ob_id].push_back(elem1);
				conB[ob_vid].push_back(elem2);
			}

			vector<unsigned>().swap(lbs);
			vector<unsigned>().swap(fls);

			// ====== 映射关系不对  =======
			vector<unsigned>& rts = root[i];

			for (int j=0; j<topk; ++j){
				// j 实际是对应 boundary graph 的id
				unsigned ob_vid = j, dis = rts[j];

				if (ob_vid == ob_id or dis == 0) continue;

				// if (oid == 5)
				// 	cout<<p2v_Global[ob_vid]<<"  *  "<<dis<<endl;

				unsigned elem1 = ob_vid << MAXMOV1 | dis,
						 elem2 =  ob_id << MAXMOV1 | dis;

				// if (_my_rank == 1)
				// 	cout<<oid<<"  "<<ob_id<<"   -   "<<p2v_Global[ob_vid]<<"  "<<ob_vid<<endl;

				conB[ob_id].push_back(elem1);
				conB[ob_vid].push_back(elem2);
			}

			vector<unsigned>().swap(rts);
		}

		delete[] label;
		delete[] flags;


		vector<vector<vector<unsigned> > > CB(_num_workers);
		
		long long Edgess, aa = 0;
		for (int i=0; i<conB.size(); ++i){
			aa += conB[i].size();
		}
		
		Edgess = all_sum_LL(aa);
		worker_barrier();

		if (_my_rank == 0)
		cout<<"V: "<<BdV<<"   E: "<<Edgess/2<<endl;

		// == 收集 CopyNum 顶点的标签信息  ==
		int start = 0, batch = 10000;
		if (start+batch > BdV)
			batch = BdV-start;

		while (1){
			for (int i=start; i<start+batch; ++i){
				
				if (conB[i].size() == 0) continue;

				conB[i].push_back(i); 
				for (int j=0; j<_num_workers; ++j)
					if (j != _my_rank)
						CB[j].push_back(conB[i]);
				conB[i].pop_back();
			}
			
			all_to_all(CB);
			worker_barrier();

			for(int i=0; i<_num_workers; ++i){
				
				for (int j=0; j<CB[i].size(); ++j){
					
					vector<unsigned>& inf = CB[i][j];
					int vid = inf[inf.size()-1];
					inf.pop_back();

					conB[vid].insert(conB[vid].end(), inf.begin(), inf.end());
					
					vector<unsigned>().swap(inf);
				}
			}
			
			for (int i=0; i<_num_workers; ++i)
				vector<vector<unsigned> >().swap(CB[i]);

			vector<vector<vector<unsigned> > >(CB).swap(CB);

			start += batch;

			int cntt = all_sum(start);

			// if (_my_rank == 0)
			// 	cout<<"Cnt: "<<cntt<<"  "<<BdV<<endl;

			if (cntt == _num_workers*BdV) break;

			if (start+batch > BdV)
				batch = BdV-start;
		}


		for (int i=0; i<conB.size(); ++i){
			
			if (conB[i].size() == 0) continue;

			vector<unsigned>& edges = conB[i];
			sort(edges.begin(), edges.end());
			int p = 1;
			for( int j = 1; j < (int) edges.size(); ++j ){
				unsigned id1 = edges[j-1] >> MAXMOV1,
						 id2 = edges[j] >> MAXMOV1;
				if( id1 != id2 ) 
					edges[p++] = edges[j];
			}
			edges.resize(p);
		}

		// ========= conB 格式转化 ===========
		MAXMOV = MAXMOV1;
		MASK = MASK1;

		vector<vector<unsigned> >().swap(con);
		con.resize(BdV);

		for (int i=0; i<conB.size(); ++i){
			
			vector<unsigned>& local = conB[i];
			
			for( int p = 0; p < local.size(); ++p ) {

				unsigned id = local[p] >> MAXMOV1, 
						 dd = local[p] & MASK1;
				
				unsigned elem = i << MAXMOV1 | dd;

				con[id].push_back(elem);
			}
		}
		
		vector<vector<unsigned> >().swap(conB);
	
		for (int i=0; i<con.size(); ++i){
			vector<unsigned>(con[i]).swap(con[i]);
		}

		worker_barrier();
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
		topk = 5;

		const char *part_path = params.partition_path.c_str(); 
		infile.open(part_path);
		string s;
		if(!infile.is_open()){
			cout<<"No such file!"<<endl;
			exit(-1);
		}

		int nnn = 0;
		ActiveV.resize(_num_workers);

		while(getline(infile, s)){
			char* part = new char[strlen(s.c_str())+1];
			strcpy(part, s.c_str());
			int pt = atoi(part);
			v2part.push_back( atoi(part) ); // 
			ActiveV[pt] += 1;
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

		BdV = all_sum(BoundNum);
		totalV = all_sum(vertexes.size());
		ALLVertices = all_sum(vertexes.size());
		long long totalE = all_sum_LL(totalEdges);

		// cout<<"rank: "<<_my_rank<<"  V: "<<vertexes.size()<<"   E: "<<LocaEdges<<endl;

		blockInit(vertexes, blocks); 

		if (_my_rank == 0){
			cout<<"Total V: "<<totalV<<"  E: "<<totalE/2<<endl;
			cout<<"Bound V: "<<BdV<<endl;
		}
		
		float t = omp_get_wtime();
		
		construct_bp_label();
		
		worker_barrier();

		Part2hop();

		// string new_filename = params.output_path + "In_"+ to_string(_my_rank);
		// DHI_Store(new_filename);

		for (int i=BoundNum; i<con.size(); ++i){
			vector<unsigned>().swap(label[i]);
			vector<unsigned>().swap(root[i]);
		}

		worker_barrier();
		
		float ttt1 = omp_get_wtime();
		
		BoundGraph();

		float ttt2 = omp_get_wtime();

		Core2hop();

		// new_filename = params.output_path + "Bd_"+ to_string(_my_rank);
		// DHBound_Store(new_filename);
		
		float ttt = omp_get_wtime()-t;
		if (_my_rank == 0){
			cout<<"time: "<<ttt<<" s"<<endl;
			cout<<"Inner time: "<<ttt1-t<<" s  "<<ttt2-ttt1<<"  s"<<endl;
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














