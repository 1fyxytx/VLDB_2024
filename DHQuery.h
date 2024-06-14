#include "../utils/communication.h"
#include "../blogel/BVertex.h"
#include "../blogel/Block.h"
#include "../blogel/BWorker.h"
#include "../blogel/BGlobal.h"
#include "../blogel/BType.h"
#include <queue>
#include <vector>
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

using namespace std;


struct STLCRVertexValue {
	vector<pair<unsigned, unsigned> > indL, indE; // 两种index
	int src, dst;

    // void set_Edge(unsigned val, int flg){

	// 	if (flg == 1)
	// 		indE.push_back(val);
	// 	else
	// 		indL.push_back(val);
    // }

	void Initial(){
		src = 1000, dst = -1000;
	}

	void empty(){
		vector<pair<unsigned, unsigned> >().swap(indL);
		vector<pair<unsigned, unsigned> >().swap(indE);
	}

    friend ibinstream & operator<<(ibinstream & m, const STLCRVertexValue & v){
		m<<v.indL;
		m<<v.indE;
		m<<v.src;
		m<<v.dst;

    	return m;
    }


    friend obinstream & operator>>(obinstream & m, STLCRVertexValue & v){
		m>>v.indL;
    	m>>v.indE;
		m>>v.src;
		m>>v.dst;

		return m;
    }
};


struct MsgInf {
	int dis;

	MsgInf(){dis = -1;}
    
	MsgInf(int d){dis = d;}

    friend ibinstream& operator<<(ibinstream& m, const MsgInf& idm)
    {
        m << idm.dis;

        return m;
    }

    friend obinstream& operator>>(obinstream& m, MsgInf& idm)
    {
        m >> idm.dis;

        return m;
    }
};


class LCRVertex : public BVertex<VertexID, STLCRVertexValue, MsgInf>{
public:
	virtual void compute(MessageContainer& messages){}

	void Inner_Compute(MessageContainer& messages){
		if (step_num() == 1){ // 添加元素，剩下的在block里面做
			if (id == src) 
				value().src = 0;
			if (id == dst) 
				value().dst = 0;
		}else{ // 第二个超步，将message中的顶点添加到send中就完成
			for (MessageIter it = messages.begin(); it != messages.end(); it++) {
				if ((*it).dis > 0 and value().src > (*it).dis)
					value().src = (*it).dis;
				
				if ((*it).dis < 0 and value().dst < (*it).dis)
					value().dst = (*it).dis;
			}
			
			if (value().src < 1000 and value().dst > -1000){
				if (khop > value().src-value().dst)
					khop = value().src-value().dst;
			}
		}

		vote_to_halt();
    }
};


class LCRBlock : public Block<char, LCRVertex, MsgInf> {
public:
	//
	virtual void compute(MessageContainer &messages, VertexContainer &vertexes){
		// == 先执行 内部的传递 ==

		if (_my_rank == v2part[src]){
		int place = vert2place[src];
		vector<pair<unsigned, unsigned> >& labs = vertexes[place]->value().indL;

		for (int i=0; i<labs.size(); ++i){
			int vid = labs[i].first, dis = labs[i].second;
			// cout<<vid<<"  "<<dis<<endl;

			int pp = vert2place[vid];
			vertexes[pp]->value().src = dis; // 更新

			vector<pair<unsigned, unsigned>>& outlabs = vertexes[pp]->value().indE;

			for (int j=0; j<outlabs.size(); ++j){
				int vid = outlabs[j].first, dis1 = outlabs[j].second;

				MsgInf elem(dis+dis1); // 传输
				vertexes[0]->send_message(vid, v2part[vid], elem);
			}
		}

		if (labs.size() == 0){ // === 边界点，直接第二步 ===
			vector<pair<unsigned, unsigned> >& outlabs = vertexes[place]->value().indE;

			for (int j=0; j<outlabs.size(); ++j){
				int vid = outlabs[j].first, dis = outlabs[j].second;

				MsgInf elem(dis); // 传输
				vertexes[0]->send_message(vid, v2part[vid], elem);
			}
		}

		}
// ===============================================

// ===============================================
		if (_my_rank == v2part[dst]){

int place = vert2place[dst];
vector<pair<unsigned, unsigned> >& labs = vertexes[place]->value().indL;
for (int i=0; i<labs.size(); ++i){
	int vid = labs[i].first, dis = labs[i].second;
	int pp = vert2place[vid];
	vertexes[pp]->value().dst = -dis; // 更新

	vector<pair<unsigned, unsigned> >& outlabs = vertexes[pp]->value().indE;

	for (int j=0; j<outlabs.size(); ++j){
		int vid = outlabs[j].first, dis1 = outlabs[j].second;

		MsgInf elem(-dis-dis1); // 传输
		vertexes[0]->send_message(vid, v2part[vid], elem);
	}
}

if (labs.size() == 0){ // === 边界点，直接第二步 ===
	vector<pair<unsigned, unsigned> >& outlabs = vertexes[place]->value().indE;

	for (int j=0; j<outlabs.size(); ++j){
		int vid = outlabs[j].first, dis = outlabs[j].second;

		MsgInf elem(-dis); // 传输
		vertexes[0]->send_message(vid, v2part[vid], elem);
	}
}


		}
	}

};


class LCRBlockWorker : public BWorker<LCRBlock>{
public:
	vector<vector<vector<unsigned> > > Lab1, Lab2;
	vector<int> D1, D2;
	vector<unsigned> *label, *outlabel;
	vector<pair<int, int> > QueryPairs;
	unsigned MAXDIS, MAXMOV, MASK;
	unsigned MAXDIS1, MAXMOV1, MASK1;

	virtual LCRVertex *toVertex(char *line){
		LCRVertex *v;
		return v;
	}


	virtual void blockInit(VertexContainer &vertexes, BlockContainer &blocks){
		
		vector<int> srcL, dstL;

		for (int i=0; i<10000; ++i){
			int vid = vertexes[i]->id;
			if (_my_rank % 2 == 1)  srcL.push_back(vid);
			else                    dstL.push_back(vid);
		}

		vector<vector<int> > SrcL(_num_workers), DstL(_num_workers);
		for (int i=0; i<_num_workers; ++i){
			SrcL[i] = srcL; DstL[i] = dstL;
		}

		all_to_all_cat(SrcL, DstL);
		worker_barrier();

		srcL.clear(), dstL.clear();

		for (int i=0; i<_num_workers; ++i){
			srcL.insert(srcL.end(), SrcL[i].begin(), SrcL[i].end());
			dstL.insert(dstL.end(), DstL[i].begin(), DstL[i].end());
		}

		vector<vector<int> >().swap(SrcL);
		vector<vector<int> >().swap(DstL);

		sort(srcL.begin(), srcL.end());
		sort(dstL.begin(), dstL.end());

		for (int i=0; i<srcL.size(); ++i){
			if (srcL[i] != dstL[i]){
				pair<int, int> elem(srcL[i], dstL[i]);
				QueryPairs.push_back(elem);
			}
		}

	}


	void active_LCR_vcompute(){
		active_vcount = 0;
		VMessageBufT* mbuf = (VMessageBufT*)get_message_buffer();
		vector<MessageContainerT>& v_msgbufs = mbuf->get_v_msg_bufs();
		for (BlockIter it = blocks.begin(); it != blocks.end(); it++)
		{
			LCRBlock* block = *it;
			for (int i = block->begin; i < block->begin + block->size; i++)
			{
				if (v_msgbufs[i].size() == 0)
				{
					if (vertexes[i]->is_active())
					{
						block->activate(); //vertex activates its block
						vertexes[i]->Inner_Compute(v_msgbufs[i]);
						if (vertexes[i]->is_active())
							active_vcount++;
					}
				}
				else
				{
					block->activate(); //vertex activates its block
					vertexes[i]->Inner_Compute(v_msgbufs[i]);
					if (vertexes[i]->is_active())
						active_vcount++;
				}
			}
		}
	}


	void all_LCR_vcompute(){
		active_vcount = 0;
		VMessageBufT* mbuf = (VMessageBufT*)get_message_buffer();
		vector<MessageContainerT>& v_msgbufs = mbuf->get_v_msg_bufs();
		for (BlockIter it = blocks.begin(); it != blocks.end(); it++)
		{
			LCRBlock* block = *it;
			block->activate(); //vertex activates its block
			for (int i = block->begin; i < block->begin + block->size; i++)
			{
				vertexes[i]->activate();
				vertexes[i]->Inner_Compute(v_msgbufs[i]);
				v_msgbufs[i].clear(); //clear used msgs
				if (vertexes[i]->is_active())
					active_vcount++;
			}
		}
	}


	void all_bcompute_LCR(){
		active_bcount = 0;
		BMessageBufT* mbuf = (BMessageBufT*)get_bmessage_buffer();
		vector<BMessageContainerT>& b_msgbufs = mbuf->get_b_msg_bufs();
		for (int i = 0; i < blocks.size(); i++){
			blocks[i]->activate();
			blocks[i]->compute(b_msgbufs[i], vertexes);
			b_msgbufs[i].clear(); //clear used msgs
			if (blocks[i]->is_active())
				active_bcount++;
		}
	};


	void InLRead(string indexName){
		label = new vector<unsigned>[totalV];
		FILE *fin = fopen( (indexName).c_str(), "rb" );
		int n;
		fread(&n, sizeof(int), 1, fin);
		int *len = new int[n];
		fread(len, sizeof(int), n, fin); // 读n条数据

		for( int i = 0; i < n; ++i ) {
			unsigned *s = new unsigned[len[i]];

			fread(s, sizeof(unsigned), len[i], fin);

			int vid = s[0];
			label[vid].reserve(len[i]);
			label[vid].assign(s, s + len[i]);

			delete[] s;
		}
	}


	void BdLRead(string indexName){
		outlabel = new vector<unsigned>[totalV];
		FILE *fin = fopen( (indexName).c_str(), "rb" );
		int n;
		fread(&n, sizeof(int), 1, fin);
		int *len = new int[n];
		fread(len, sizeof(int), n, fin); // 读n条数据

		for( int i = 0; i < n; ++i ) {
			unsigned *s = new unsigned[len[i]];

			fread(s, sizeof(unsigned), len[i], fin);

			unsigned vid = s[len[i]-2];
			
			if (vid > totalV or vid <= 0) 
				cout<<vid<<"  "<<totalV<<endl;

			outlabel[vid].reserve(len[i]-2);
			outlabel[vid].assign(s, s + len[i]-2);

			delete[] s;
		}
	}


	void LabelReSet(){
		for (int i=0; i<totalV; ++i){
			if (v2part[i] == _my_rank){
				LCRVertex* v = new LCRVertex;
				v->id = i;
				v->value().src = 1000;
				v->value().dst = -1000;
				v->bid = 0;
				load_vertex(v);
				vert2place[i] = vertexes.size()-1;
			}

			if (label[i].size() > 0){ // interior vertex
				for(int kk=0; kk<(int)(label[i].size()/2); ++kk){
					unsigned vvid = label[i][2*kk], ddis = label[i][2*kk+1];
					pair<unsigned, unsigned> elem(vvid, ddis);
					vertexes[vert2place[i]]->value().indL.push_back(elem);
				}
			}else{ // == clear redundant labels ==
				int p1 = 0;
				for (int j=0; j<int(outlabel[i].size()/2); ++j){
					unsigned vid = outlabel[i][2*j];
					if (vid % _num_workers == _my_rank){
						outlabel[i][p1++] = outlabel[i][2*j];
						outlabel[i][p1++] = outlabel[i][2*j+1];
					}
				}
				outlabel[i].resize(p1);
				outlabel[i].push_back(i);
				Lab1[v2part[i]].push_back(outlabel[i]); 
			}
		}	

		all_to_all(Lab1);
		worker_barrier();

		for (int i=0; i<_num_workers; ++i){
			for (int j=0; j<Lab1[i].size(); ++j){
				vector<unsigned>& lab = Lab1[i][j];
				int vid = lab[lab.size()-1];
				lab.pop_back();	

				vector<pair<unsigned, unsigned> >& edge = vertexes[vert2place[vid]]->value().indE;
				for(int kk=0; kk<(int)(lab.size()/2); ++kk){
					unsigned vvid = lab[2*kk], ddis = lab[2*kk+1];
					pair<unsigned, unsigned> elem(vvid, ddis);
					edge.push_back(elem);
				}
			}
		}

		vector<vector<vector<unsigned> > >().swap(Lab1);
		delete[] label;
		delete[] outlabel;
	}


	void run_LCR(const WorkerParams& params){
		ifstream infile, infile1, infile2; // 先把分区文件读取
		
		src = params.src, dst = params.dst;
		khop = 10000;

		const char *part_path = params.partition_path.c_str(); 
		infile.open(part_path);
		string ss;
		if(!infile.is_open()){
			cout<<"No such file!"<<endl;
			exit(-1);
		}

		int nnn = 0;
		vector<int> ActiveV(_num_workers);

		while(getline(infile, ss)){
			char* part = new char[strlen(ss.c_str())+1];
			strcpy(part, ss.c_str());
			v2part.push_back( atoi(part) );
			ActiveV[atoi(part)] += 1;
			v2degree.push_back(0);
			delete part;
			nnn += 1;
		}

		totalV = nnn;
		Lab1.resize(_num_workers), Lab2.resize(_num_workers);
        
		string indexName = params.input_path+"In_"+ to_string(_my_rank);
		InLRead(indexName);

		indexName = params.input_path+"Bd_"+ to_string(_my_rank);
		// if (_my_rank == 0)
		BdLRead(indexName);

		worker_barrier();

// 	    // ==== label and outlabel 的元素转移到顶点中  ====
		LabelReSet();

		blockInit(vertexes, blocks);

		int prev = -1;
		LCRBlock* block = NULL;
		int pos;
		for (pos = 0; pos < vertexes.size(); pos++){
			int bid = vertexes[pos]->bid; 
			if (bid != prev){
				if (block != NULL){
					block->size = pos - block->begin;
					blocks.push_back(block);
				}
				block = new LCRBlock;
				prev = block->bid = bid;
				block->begin = pos;
			}
		}

        if (block != NULL){
            block->size = pos - block->begin;
            blocks.push_back(block);
        }
        // active_bcount = getBNum(); //initially, all blocks are active


        get_bnum() = all_sum(getBNum());
        get_vnum() = all_sum(getVNum());

        // blockInit(vertexes, blocks); //setting user-defined block fields, 用户可以指定给block构建什么信息

        vmessage_buffer->init(vertexes);
        bmessage_buffer->init(blocks);

        worker_barrier();

// =================================================================
		float TotalTime = 0;
		long long TotalMessage = 0;
		unordered_map<int, float> time2d, comm2d;
		int cnt1 = 0;
		
        if (compute_mode == VB_COMP){
		for (int ii=0; ii<QueryPairs.size(); ++ii){

			long long step_vmsg_num;
			long long step_bmsg_num;
			long long global_vmsg_num = 0;
			long long global_bmsg_num = 0;

			src = QueryPairs[ii].first, dst = QueryPairs[ii].second;

        	float s_1 = clock();

        	global_step_num = 0;
        	while (true){
				global_step_num++;

				if (global_step_num == 1){
					all_LCR_vcompute();
				}else{
					active_LCR_vcompute();
					break;
				}

				all_bcompute_LCR();

				// ---  统计传输数据量  -----
                step_vmsg_num = master_sum_LL(vmessage_buffer->get_total_msg());
                step_bmsg_num = master_sum_LL(bmessage_buffer->get_total_msg());
                if (_my_rank == 0) {
                    global_vmsg_num += step_vmsg_num;
                    global_bmsg_num += step_bmsg_num;
                }

                vmessage_buffer->sync_messages();
                bmessage_buffer->sync_messages();

                //===================
                worker_barrier();
        	}

			worker_barrier();
        	float s_2 = clock();

			TotalTime += (s_2-s_1);
			TotalMessage += global_vmsg_num;

			if (_my_rank == 0 and ii % 20 == 0){
			cout<<" Cnt: "<<ii<<"  Average: "<<float(TotalTime)/(ii*CLOCKS_PER_SEC)<<" s , Comm: "<<float(TotalMessage)/(ii*1024)<<endl;

			}

			khop = 10000;
			for (int i=0; i<vertexes.size(); ++i)
				vertexes[i]->value().Initial();
		}

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

