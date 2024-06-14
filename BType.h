#ifndef BTYPE_H_
#define BTYPE_H_
#include "../utils/type.h"
#include <bitset>
#include <queue>
#define MAXD 120
#define MAXINT ((unsigned) 4294967295)
#define N_ROOTS 0

struct triplet {
    VertexID vid;
    int bid;
    int wid;

    bool operator==(const triplet& o) const
    {
        return vid == o.vid;
    }

    friend ibinstream& operator<<(ibinstream& m, const triplet& idm)
    {
        m << idm.vid;
        m << idm.bid;
        m << idm.wid;
        return m;
    }
    friend obinstream& operator>>(obinstream& m, triplet& idm)
    {
        m >> idm.vid;
        m >> idm.bid;
        m >> idm.wid;
        return m;
    }
};

struct Index_Element{
	int id;
	int labels;

	Index_Element(int idl, int label){
		id = idl;
		labels = label;
	}
};

struct Bound_Element{
	// 针对单边多标签的组成
	int id;
	int labels;
	int adj_label;

	Bound_Element(int id1, int label, int adjL){
		id = id1;
		labels = label;
		adj_label = adjL;
	}
};

struct DSR_Element{
	int id;
	set<int> infs;

	DSR_Element(int id1){
		id = id1;
	}

	DSR_Element(int id1, set<int> inf1){
		id = id1;
		infs = inf1;
	}

	DSR_Element(int id1, int inf1){
		id = id1;
		infs.insert(inf1);
	}
};

int BitCount5(unsigned int n){
	int ii = 1;
	if (ii==1){
		unsigned int tmp = n - ((n >> 1) & 033333333333) - ((n >> 2) & 011111111111);
		return ((tmp + (tmp >> 3)) & 030707070707) % 63;
	}else
		return n;
}

struct Index{
	int labels;

	Index(){
	}

	Index(int label){
		labels = label;
	}

	bool operator < (const Index&a) const{
		if (BitCount5(labels) < BitCount5(a.labels)){
			return true;
		}
		else if (BitCount5(labels) == BitCount5(a.labels)){
			if (labels < a.labels){
				return true;
			}
			return false;
		}else
			return false;
	}

	bool operator > (const Index&a) const{
		if (BitCount5(labels) > BitCount5(a.labels)){
			return true;
		}
		else if (BitCount5(labels) == BitCount5(a.labels)){
			if (labels > a.labels){
				return true;
			}
			return false;
		}else
			return false;
	}

	void operator = (const Index&a){
		labels = a.labels;
	}

	bool operator == (const Index&a) const{
		if (labels == a.labels)
			return true;
		else
			return false;
	}

	friend ibinstream& operator<<(ibinstream& m, const Index& idm)
	{
		m << idm.labels;

		return m;
	}

	friend obinstream& operator>>(obinstream& m, Index& idm)
	{
		m >> idm.labels;

		return m;
	}

};

struct Dynamic_Element{
	// 针对单边多标签的组成
	int id; // 激活点
	int labels; // 当前标签
	int adj_label; // 邻接边的标签
	int flag;  // 确认对应邻边，1-新边；0-所有边

	Dynamic_Element(int id1, int label, int adjL, int ff){
		id = id1;
		labels = label;
		adj_label = adjL;
		flag = ff;
	}
};

int totalLength = 0;
int totalQuery = 0;
long totalEdges = 0;
int totalV = 0, Nums = 0, ALLVertices=0;
int edgeNum = 0, labelnum = 6;
int cnt, boundV = 0, BdV = 0;
float time11, time2;
vector<double> time12;
vector<int> In_Boundary, Out_Boundary, srcList, dstList;

double totalSize = 0;

int convert_to_binary_label(int temp)
{
	if (temp == 0 || temp == 1)
	{
		return temp+1;
	}
	int result = 1;
	result = result << (temp);
	return result;
}


//int VerdicateSubset(vector<int> label_vector, int label){ // 是子集，返回1
//	int res = 0;
//
//	for (int i=0; i<label_vector.size(); ++i){
//		cnt += 1;
//		if (label_vector[i] > label)
//			break;
//
//		int r = label_vector[i]&label;
//		if ( r == label_vector[i]){
//			return 1;
//		}
//	}
//
//	return res;
//}


int VerdicateSubset(vector<Index>& label_vector, int label){ // 是子集，返回1
	int res = 0;

	for (int i=0; i<label_vector.size(); ++i){
		cnt += 1;

		int r = label_vector[i].labels & label;
		if ( r == label_vector[i].labels){
			return 1;
		}

		if (BitCount5(label_vector[i].labels) > BitCount5(label))
			break;
	}

	return res;
}

int lc = 15, numbers = 10000;
vector<vector<int> > Src, Dst, src1, dst1;

int Query_Result(map<int, vector<Index> >& src, int active_vid,
		map<int, vector<Index> >& dst, int target_vid, int label){
	// src 的排序值肯定在 dst 之上

	int res = 1, res1=0, res2=0;
	int mid_v;
	map<int, vector<Index> >::iterator it;

	if ( src.size()==0 && dst.size()==0 )
		return 1;
	else{
		if (dst.size()>0){
			if (dst.find(active_vid) == dst.end()){

			}
			else{ // dst[active_vid]
			    res = VerdicateSubset(dst[active_vid], label);
				if (res == 1) // 新标签不需要更新
					return 0;
				else
					res = 1;
			}
		}// src不会有target_vid的直连标签

		if ( src.size()==0 || dst.size()==0 )
			return 1;

		if ( src.size() < dst.size() )// 从 src中的标签点进行遍历，判定是否满足查询程序
			for ( it=src.begin(); it!=src.end(); ++it ){
				if ( dst.find(it->first) == dst.end() )
					continue;

				mid_v = it->first;
				// 如果同时在 src[middle_v]和dst[middle_v]两者中都找到了label的子集，那么label这个标签就不需要添加到索引中
				res1 = VerdicateSubset(src[mid_v], label);
				if (res1 == 1) // 表示在src这边找到了现有标签的子集，需要继续找dst的子集
					res2 = VerdicateSubset(dst[mid_v], label);

				if (res2 == 1)
					return 0; // 现有标签满足，不用更新
			}
		else
			for (it=dst.begin(); it!=dst.end(); ++it){
				if ( src.find(it->first) == src.end() )
					continue;

				mid_v = it->first;
				res1 = VerdicateSubset(dst[mid_v], label);
				if (res1 == 1) // 表示在src这边找到了现有标签的子集，需要继续找dst的子集
					res2 = VerdicateSubset(src[mid_v], label);

				if (res2 == 1)
					return 0; // 现有标签满足，不用更新
			}
	}
	// 没有在中间返回结果，那么表示这个标签需要进行更新

//	float ss_2 = clock();
//	time12 += ss_2 - ss_1;

	return res;
}

int Query_Result_2(vector<int>& L1, vector<int>& L2){
	
	int i=0, j=0, end1=L1.size(), end2=L2.size();

	if (end1==0 || end2==0)
		return 1;

	if (L1[0]>L2[end2-1] || L1[end1-1]<L2[0])
		return 1;

	while(i<end1 && j<end2){
		if (L1[i]==L2[j])
			return 0;
		else if (L1[i] < L2[j])
			i++;
		else
			j++;
	}

	return 1;
}

vector<Index> InnerPathVerdicate(map<int, vector<Index> >& src, int src_vid,
		map<int, vector<Index> >& dst, int dst_vid){

	map<int, vector<Index> >::iterator it;
	set<Index>::iterator its;
	set<Index> LabelSet;
	vector<Index> ress;
	ress.resize(0);

	for ( it=src.begin(); it!=src.end(); ++it ){
		if ( dst.find(it->first) == dst.end() )
			continue;

		vector<Index>& L1 = src[it->first];
		vector<Index>& L2 = dst[it->first];

		for (int ii=0; ii<L1.size(); ++ii){
			for (int jj=0; jj<L2.size(); ++jj){
				Index elem(L1[ii].labels | L2[jj].labels);
				LabelSet.insert(elem);
			}
		}
	}

	if (src.find(dst_vid) != src.end()){
		LabelSet.insert(src[dst_vid].begin(), src[dst_vid].end());
	}

	if (dst.find(src_vid) != dst.end()){
		LabelSet.insert(dst[src_vid].begin(), dst[src_vid].end());
	}

	// 已经完成了去重，下一步就是要去冗余边
	for (its=LabelSet.begin(); its!=LabelSet.end(); ++its){

		int res = VerdicateSubset(ress, (*its).labels);

		if (res == 1){
			Index newOne((*its).labels);
			ress.insert( lower_bound(ress.begin(), ress.end(), newOne), newOne );
		}

	}

	return ress;
}


set<int> convert_to_label_set(int temp)
{
	set<int> result;
	int i = 0;
	while (temp != 0)
	{
		if ((temp & 1) == 1)
		{
			result.insert(convert_to_binary_label(i));
		}
		i++;
		temp = temp >> 1;
	}
	return result;
}

vector<int> L1, R1, F1;

vector<int> v2part, v2degree, v2bound, vFlag, v2degree_R; // 全局描述变量

vector<pair<int, int> > newEdges;



int global_flag = 1;
float query_time = 0, insert_time = 0, label_time = 0, ss1 = 0, ss2=0;
set<int> labelNum;

set<int> Node_B, Node_I2B, Node_I; // 所有从inner变成boundary的顶点
vector<set<int> > vertices_B, vertices_I2B;
vector<pair<int, int> > sortList, sortTotal, sortTotal_R;
map<int, queue<int> > Update_In_Oneplus, Update_Out_Oneplus, Bound_In_Oneplus, Bound_Out_Oneplus; // 存储标签数量+1的队列元素


// 第二阶段的程序，暂时不改

void PartRevise(vector<int>& part, int num){
	int totalV = part.size();
	int maxV = int(totalV/num) + 1;

	for (int i=0; i<totalV; i++){
		part[i] = i / maxV;
	}
}

#endif

// =================================================
// ======== 4-th paper ===  temporal graph  =======

struct tmptriple{
	int id;
	int st;
	int et;

	tmptriple(){
		st=-1, et=-1, id=-1;
	}

	tmptriple(int vid, int stime, int etime){
		st = stime; et = etime; id = vid;
	}

	bool operator < (const tmptriple&a) const{
		if (st > a.st)
			return true;
		else if (st == a.st and (et-st) < (a.et-a.st))
			return true;
		else if (st == a.st and (et-st) == (a.et-a.st) and id<a.id)
			return true;
		else
			return false;
	}

	bool operator > (const tmptriple&a) const{
		if (st < a.st)
			return true;
		else if (st == a.st and (et-st) > (a.et-a.st))
				return true;
		else if (st == a.st and (et-st) == (a.et-a.st) and id > a.id)
			return true;
		else
			return false;
	}

	void operator = (const tmptriple&a){
		st = a.st; et = a.et; id = a.id;
	}


	friend ibinstream& operator<<(ibinstream& m, const tmptriple& idm)
	{
		m << idm.st; m << idm.et; m << idm.id;

		return m;
	}

	friend obinstream& operator>>(obinstream& m, tmptriple& idm)
	{
		m >> idm.st; m >> idm.et; m >> idm.id;

		return m;
	}
};


struct tmp{
	int st;
	int et;

	tmp(){
		st=-1, et=-1;
	}

	tmp(int stime, int etime){
		st = stime; et = etime;
	}

	bool operator < (const tmp&a) const{
		if (st > a.st)
			return true;
		else if (st == a.st and (et-st) < (a.et-a.st))
			return true;
		else
			return false;
	}

	bool operator > (const tmp&a) const{
		if (st < a.st)
			return true;
		else if (st == a.st and (et-st) > (a.et-a.st))
			return true;
		else
			return false;
	}

	bool operator == (const tmp&a) const{
		if (st == a.st and (et-st) == (a.et-a.st))
			return true;
		else
			return false;
	}

	void operator = (const tmp&a){
		st = a.st; et = a.et;
	}

	friend ibinstream& operator<<(ibinstream& m, const tmp& idm)
	{
		m << idm.st; m << idm.et;

		return m;
	}

	friend obinstream& operator>>(obinstream& m, tmp& idm)
	{
		m >> idm.st; m >> idm.et;

		return m;
	}
};


struct tmplus{
	tmp interval;
	vector<int> idList;

	tmplus(){
		interval.st=-1, interval.et=-1; idList.resize(0);
	}

	tmplus(int stime, int etime, vector<int>& ids){
		interval.st = stime; interval.et = etime; idList=ids;
	}

	bool operator < (const tmplus& a) const{
		if (interval < a.interval)
			return true;
		else // >=的情况合并一起
			return false;
	}

	bool operator > (const tmplus& a) const{
		if (interval > a.interval)
			return true;
		else // >=的情况合并一起
			return false;
	}

	bool operator == (const tmplus& a) const{
		if (interval == a.interval)
			return true;
		else // >=的情况合并一起
			return false;
	}
};


struct tmpOut{
	int st;
	int et;

	tmpOut(){
		st = -1; et = -1;
	}

	tmpOut(int stime, int etime){
		st = stime; et = etime;
	}

	bool operator < (const tmpOut&a) const{
		if (st < a.st)
			return true;
		else if (st == a.st and (et-st) < (a.et-a.st))
			return true;
		else
			return false;
	}

	bool operator > (const tmpOut&a) const{
		if (st > a.st)
			return true;
		else if (st == a.st and (et-st) > (a.et-a.st))
				return true;
		else
			return false;
	}

	bool operator == (const tmpOut&a) const{
		if (st == a.st and (et-st) == (a.et-a.st))
				return true;
		else
			return false;
	}

	void operator = (const tmpOut&a){
		st = a.st; et = a.et;
	}


	friend ibinstream& operator<<(ibinstream& m, const tmpOut& idm)
	{
		m << idm.st; m << idm.et;
		return m;
	}

	friend obinstream& operator>>(obinstream& m, tmpOut& idm)
	{
		m >> idm.st; m >> idm.et;
		return m;
	}
};


struct tmpIn{
	int st;
	int et;

	tmpIn(){
		st = -1; et = -1;
	}

	tmpIn(int stime, int etime){
		st = stime; et = etime;
	}

	bool operator < (const tmpIn&a) const{
		if (et > a.et)
			return true;
		else if (et == a.et and (et-st) < (a.et-a.st))
			return true;
		else
			return false;
	}

	bool operator > (const tmpIn&a) const{
		if (et < a.et)
			return true;
		else if (et == a.et and (et-st) > (a.et-a.st))
				return true;
		else
			return false;
	}

	bool operator == (const tmpIn&a) const{
		if (et == a.et and (et-st) == (a.et-a.st))
			return true;
		else
			return false;
	}

	void operator = (const tmpIn&a){
		st = a.st; et = a.et;
	}


	friend ibinstream& operator<<(ibinstream& m, const tmpIn& idm)
	{
		m << idm.st; m << idm.et;
		return m;
	}

	friend obinstream& operator>>(obinstream& m, tmpIn& idm)
	{
		m >> idm.st; m >> idm.et;
		return m;
	}
};


int Dominate(vector<tmpOut>& tmpOutV, vector<tmpIn>& tmpInV, pair<int, int> tvalue, int flag){ // 是子集，返回1
	// 当存在直连标签时，执行三种操作
	// 1.新标签可以被已有标签dominate，则返回1
	// 2.新标签可以dominate现有标签，要获取所有需要被dominated的元素，然后统一删除。
	// 3.新标签和现有标签之间无1,2关系，只返回0，还需要进一步和中间标签进行判定.
	int res = 0;
	int st = tvalue.first, et = tvalue.second;

	if (flag == 0){ // 验证tmpInV是否已经包含 tvalue
		queue<tmpIn> Del_Elems;

		for(int i=0; i<tmpInV.size(); ++i){
			int stExist = tmpInV[i].st, etExist = tmpInV[i].et;

			if (etExist <= et and stExist >= st) return 1; // 表示情况1

			if ( (etExist>et and stExist<=st) or
				 (etExist>=et and stExist<st) ) Del_Elems.push(tmpInV[i]); // 情况2，记录下需要被删除的值.

			if (etExist <= st) break; // 基于tmpIn的一个剪枝方法，此时可以直接得到情况3，无需接着遍历
		}

		if (Del_Elems.size() > 0){ // 先剔除要删除的元素，然后返回-1
			tmpIn nowOne = Del_Elems.front();  Del_Elems.pop();

			for(vector<tmpIn>::iterator iter=tmpInV.begin(); iter!=tmpInV.end();  ++iter){
				 if( *iter == nowOne){
					iter = tmpInV.erase(iter);  iter -- ;
					if (Del_Elems.size()==0)  break;
					nowOne = Del_Elems.front();  Del_Elems.pop();
				 }
			}
			return -1;
		}

	}else{
		queue<tmpOut> Del_Elems;

		for(int i=0; i<tmpOutV.size(); ++i){
			int stExist = tmpOutV[i].st, etExist = tmpOutV[i].et;

			if (etExist <= et and stExist >= st) return 1; // 表示情况1

			if ( (etExist>et and stExist<=st) or
				 (etExist>=et and stExist<st) ) Del_Elems.push( tmpOutV[i] ); // 情况2

			if (stExist >= et) break; // 情况3
		}

		if (Del_Elems.size() > 0){ // 先剔除要删除的元素，然后返回-1

			tmpOut nowOne = Del_Elems.front();  Del_Elems.pop();

			for(vector<tmpOut>::iterator iter=tmpOutV.begin(); iter!=tmpOutV.end();  ++iter){
			     if( *iter == nowOne){
			        iter = tmpOutV.erase(iter);  iter -- ;
			        if (Del_Elems.size()==0)     break;
			        nowOne = Del_Elems.front();  Del_Elems.pop();
			     }
			}

			return -1;
		}

	}

	return res;
}


int Verdicate_temporal(vector<tmpOut>& tmpOutV, vector<tmpIn>& tmpInV, pair<int, int> tvalue){
	int st = tvalue.first, et = tvalue.second;

	if (tmpOutV[tmpOutV.size()-1].st<st or tmpInV[tmpInV.size()-1].et>et)
		return 0; // 当前中间点的构建索引没有满足dominate新标签的能力

	// == 这步开始，表明至少两方都有满足约束的中间索引  ===
	for (int i=tmpOutV.size()-1; i>-1; --i){
		// 新标签肯定不会被dominated，这里就不用管是否能dominate现有的
		if (tmpOutV[i].st < st) break;

		int eti = tmpOutV[i].et;

		for (int j=tmpInV.size(); j>-1; --j){
			if (tmpInV[j].et > et) break; // 在这之后已经不需要查了。

			int stj = tmpInV[j].st;

			if (eti <= stj) // 当前情况，表明新的label可以被现有的dominate，所以可以直接返回1
				return 1;
		}
	}

	return 0;
}



int Query_Temporal(map<int, vector<tmpOut> >& Lout, int active_vid,
		map<int, vector<tmpIn> >& Lin, int flag, pair<int, int>& TLabel){
	// == 返回1，表示需要添加；返回0，表示不需要添加  ===
	int res, mid_v;
	vector<tmpOut> outL;
	vector<tmpIn> inL;

	if ( Lout.size()==0 and Lin.size()==0 )  return 1;
	else{ // flag=0表示这是一次正向遍历，active_vid针对的是Lin
		if (Lin.size() > 0 and flag == 0){

			if (Lin.find(active_vid) == Lin.end()){ }
			else{
				res = Dominate(outL, Lin[active_vid], TLabel, 0);
				if (res == 1) return 0;// 新标签不需要更新
				if (res == -1) return 1; // 需要添加新标签，已经完成内部多余元素的清除
				// 直连标签无法判断，还需要借由中间连接标签进行判断
			}
		}
		// flag=1表示这是一次反向遍历，active_vid针对的是Lout
		if (Lout.size() > 0 and flag == 1){
			if ( Lout.find(active_vid)==Lout.end() ){}
			else{
				res = Dominate(Lout[active_vid], inL, TLabel, 1);
				if (res == 1)  return 0; // 新标签不需要更新
				if (res == -1) return 1; // 需要添加新标签，同时还要在内部进行更新
				// 直连标签无法判断，还需要借由中间连接标签进行判断
			}
		}

		if (Lout.size()==0 or Lin.size()==0) return 1;
			// 经历过直连边的判定，现在需要用中间索引判断，此时若任何一方为0，则中间点无需继续判定

		if ( Lout.size() < Lin.size() ){// 从 src中的标签点进行遍历，判定是否满足查询程序
			map<int, vector<tmpOut> >::iterator it;
			for ( it=Lout.begin(); it!=Lout.end(); ++it ){
				if ( Lin.find(it->first) == Lin.end() )
					continue;

				mid_v = it->first;  //	基于 Lout[mid_v]和Lin[mid_v]进行判断。
				res = Verdicate_temporal(Lout[mid_v], Lin[mid_v], TLabel);

				if (res == 1) return 0; // 返回1表示，新的index没有必要更新到Lin或Lout中
				// 返回值为0时，需要继续遍历
			}
		} else{
			map<int, vector<tmpIn> >::iterator it;
			for (it=Lin.begin(); it!=Lin.end(); ++it){
				if ( Lout.find(it->first) == Lout.end() )
					continue;

				mid_v = it->first;
				// Verdicate_temporal返回0，表示需要继续判断；返回1，表示可以中断了
				res = Verdicate_temporal(Lout[mid_v], Lin[mid_v], TLabel);

				if (res == 1) return 0;// 返回1表示，新的index没有必要更新到Lin或Lout中
			}
		}
	}
	// 没有在中间返回结果，那么表示这个标签需要进行更新
	return 1;
}


// =========== 路径枚举的相关参数 =============
unsigned MAXDIS, MAXMOV, MASK;

void initial(long totalV){
	MAXDIS = 2; MAXMOV = 1;
	while( MAXINT / (totalV * 2) >= MAXDIS ){
		MAXDIS *= 2;
		++MAXMOV;
	}

	MASK = MAXDIS - 1;
}

int Query_PathEmu(vector<unsigned>& L1, vector<unsigned>& L2, int distance){
	int i=0, j=0, end1=L1.size(), end2=L2.size();

	while(i<end1 && j<end2){
		int id_i = L1[i] >> MAXMOV, id_j = L2[j] >> MAXMOV;

		if (id_i==id_j){
			// == 有相同邻居点时，判断distance和合成的距离的大小 ==
			int d1 = L1[i] & MASK, d2 = L2[j] & MASK;
			int sumV = d1 + d2;
//			cout<<"sumV:  "<<sumV<<endl;
			if (distance >= sumV)
				return 0;
			else
				i++, j++;
		}
		else if (id_i < id_j)
			i++;
		else
			j++;
	}

	return 1;
}

// ============== DisOracle =============
map<int, int> v2NewId, new2place;
vector<vector<int> > BoundV;
int src, dst, khop;
set<int> boundSet;
map<unsigned, unsigned> srcdis, dstdis;
vector<unsigned> src_dis, dst_dis; // 记录src(dst)和边界点之间的距离，
vector<int> v2p, p2v, v2p_Global, p2v_Global;
long long paths=0;
vector<pair<unsigned, unsigned> > distance_;
int k1, k2;
int srcmin = 999, dstmin = 999;
int maxhop=0;
int active_nums = 0, active_Edges = 0, active_tasks = 0;

// ======== Distributed enumeration 0428  =========
struct PTnode{
	int id;    // 瀵瑰簲鍥句腑鐨刬d
	// OPEN 1, CLOSED 2, POTENTIAL 3, POSTPONE 4, NULL 5
	int state; // 鍏辨湁5涓暟
	unordered_map<int, int > pt, dis; // pruning thresholds for vertices
	int parent; // 鍦ㄨ矾寰勪腑鐨勪綅缃
	int child; // 瀛╁瓙鑺傜偣鏁伴噺
	int level; // 鍦ㄦ爲涓殑灞

	PTnode(){
		id = -1, state = -1; parent = -1; child = -1; level = -1;
	}

	PTnode(int id1, int sta, int pla, int chd, int lv){
		id = id1, state = sta; parent = pla; child = chd; level = lv;
	}
};

vector<pair<int, int> > dist, newDis; // the distance to s and t
vector<int> flag;
vector<vector<pair<int, int> > > infVector;
vector<vector<vector<int> > > subpaths, Locals, Locals_R;
vector<vector<int> > trans_paths, transmit_paths, branches, graphs, graphs_R, gLocal, gEdge;
vector<vector<long long> > v2space, v2space_Reverse;
long long theta = 0;

vector< long long > task_load;
vector<pair<long long, vector<int> > > Divide_branches;
long long Total_Load = 0, totals;
int threads = 1;
int Lmax=0, Rmax=0;
vector<vector<int> > MFlag_Left, MFlag_Right;
vector<vector<vector<int> > > Lpaths, Rpaths; // 1-th 对应的hop数量
// 选择中间点时，需要确立传播方向，0表示从src开始传输
int md_flag = 0;
vector<int> midVertex;
vector<pair<int, int> > mid2part; // 划分到其他分区的中间点的集合
vector<long long> pth_L, pth_R;



unordered_map<int, vector<vector<int> > > P_thread, P_R_thread;
int maxlevel = 0;

vector<int> pathLen;




// ==== Enumeration 02-21 ====
vector<int> v2copy, srcVertex, tgtVertex;
vector<vector<int> > srcVList, taskQueue, transQueue;
unordered_map<int, vector<vector<int> > > MiddlePaths, Middle_Tgt_Paths; // 多加一个维度，用于hop
int kk1, kk2, kmax;
vector<vector<int> > SubP, SubPTgt;
vector<pair<long long, int > > Tasks;
vector<pair<long long, vector<int> > > Tasks1;

string new_filename1;
long long cnt4=0, cnt5=0, cnt6=0, cnt7=0;

int FLGs = 0;
unordered_map<int, vector<vector<int> > > ML_Paths, MR_Paths;
vector<int> v2Middle;
long long Wavg = 0, W1 = 0;

























// ========== Span reachability  ============
struct sEdge{
    int adj, tlb;

    sEdge(){ adj = 0, tlb = 0; }

    sEdge(int id, int lb){ adj = id, tlb = lb; }
};

struct sInf{
	int vid, tbeg, tend;

	sInf(){ vid=0, tbeg=0, tend=0; }

    sInf(int id, int t1, int t2){
    	vid = id, tbeg = t1, tend = t2;
    }
};


struct sLabel{ // index
    int tbeg, tend;

    sLabel(){ tbeg=0, tend=0; }

    sLabel( int t1, int t2){
    	tbeg = t1, tend = t2;
    }

	// == 排序方式 ==
	bool operator < (const sLabel& a) const{
		if (tend-tbeg < a.tend-a.tbeg)
			return true;
		else if(tend-tbeg == a.tend-a.tbeg){
			if (tbeg < a.tbeg) 
				return true;
			else
				return false;
		}else
			return false;
	}

	bool operator == (const sLabel& a) const{
		if (tbeg == a.tbeg and tend == a.tend)
			return true;
		else
			return false;
	}

	void operator = (const sLabel&a){
		tbeg = a.tbeg;
		tend = a.tend;
	}
};

vector<pair<int, int> > ALL_sort;

int VerdicateSpan(vector<sLabel>& V1, vector<sLabel>& V2, sLabel& tlb){

}

int Query_Span(unordered_map<int, std::vector<sLabel> >& L1, 
               unordered_map<int, std::vector<sLabel> >& L2, 
			   sLabel& tlb){
	// 检查tlb是否可以被现有的标签所表征
	for (auto it=L1.begin(); it!=L1.end(); ++it){
		if ( L2.find(it->first) == L2.end() )
			continue;

	}

}



// ===== Path counting =====
struct treeinf{ // index
    int id, dLocal;

    treeinf(){ id=0, dLocal=0; }

    treeinf( int id1, int d1){
    	id = id1, dLocal = d1;
    }

	// == 排序方式 ==
	bool operator < (const treeinf& a) const{
		if (v2degree[id] < v2degree[a.id])
			return true;
		else if(v2degree[a.id] == v2degree[id]){
			if (id < a.id) 
				return true;
			else
				return false;
		}else
			return false;
	}

	void operator = (const treeinf&a){
		id = a.id;
		dLocal = a.dLocal;
	}

	bool operator == (const treeinf&a){
		if (id == a.id and dLocal == a.dLocal)
			return true;
	}

	bool operator != (const treeinf&a){
		if (id != a.id)
			return true;
	}
};

ibinstream& operator<<(ibinstream& m, const treeinf& v)
{
	m << v.id;
	m << v.dLocal;
	return m;
}

obinstream& operator>>(obinstream& m, treeinf& v)
{
	m >> v.id;
	m >> v.dLocal;
	return m;
}

int CutEdge = 0, CopyNum = 0;

struct TreeNode{
	int parent;
	vector<treeinf> ancester;

	TreeNode(){parent=-1, ancester.resize(0);}

	TreeNode(vector<treeinf>& edge){
		parent = edge[0].id;
		ancester.insert(ancester.end(), 
						edge.begin(), 
						edge.end());
	}
};

vector<TreeNode> CoreTree;

vector<vector<pair<int, int> > > DisL;
int MAX_BP_THREADS = 8;

struct BPLabel {
    uint8_t bpspt_d[N_ROOTS];
    uint64_t bpspt_s[N_ROOTS][2];
};


struct TriElem{
	int dis, vid, pts;
	
	TriElem(){ dis = 0, vid = 0, pts = 0; }

	TriElem(int d, int v, int p){
		dis = d, vid = v, pts = p;
	}

	bool operator < (const TriElem& a) const{
		if (dis < a.dis)
			return true;
		else if(dis == a.dis){
			if (vid < a.vid) 
				return true;
			else
				return false;
		}else
			return false;
	}
};


struct CoreIndex{
	unsigned val; // id+dis
	int pts; 

	CoreIndex(){ val = 0, pts = 0; }

	CoreIndex(unsigned vv, int p){
		val = vv, pts = p;
	}

	bool operator < (const CoreIndex& a) const{
		if (val < a.val)
			return true;
		else if(val == a.val){
			if (pts < a.pts) 
				return true;
			else
				return false;
		}else
			return false;
	}

	void operator = (const CoreIndex& a){
		val = a.val, pts = a.pts;
	}
};


ibinstream& operator<<(ibinstream& m, const CoreIndex& v)
{
	m << v.val;
	m << v.pts;

	return m;
}


obinstream& operator>>(obinstream& m, CoreIndex& v)
{
	m >> v.val;
	m >> v.pts;
	
	return m;
}


struct CoreElem{
	int tgt, src, dis;

	CoreElem(){tgt=-1, src=-1, dis=-1;}

	CoreElem(int t, int s, int d){
		tgt = t, src = s, dis = d;
	}

	bool operator < (const CoreElem& a) const{
		if (tgt < a.tgt)
			return true;
		else if(tgt == a.tgt){
			if (src < a.src) 
				return true;
			else if (src == a.src){
				if (dis < a.dis)
					return true;
				else
					return false;
			}
		}else
			return false;
	}
};


ibinstream& operator<<(ibinstream& m, const CoreElem& v)
{
	m << v.tgt;
	m << v.src;
	m << v.dis;

	return m;
}


obinstream& operator>>(obinstream& m, CoreElem& v)
{
	m >> v.tgt;
	m >> v.src;
	m >> v.dis;
	
	return m;
}

