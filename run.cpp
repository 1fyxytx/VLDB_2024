#include "Distri4hop.h"

using namespace std;

int main(int argc, char* argv[]){
	init_workers();
	SetLCR_Construction(  "/home/hnu/Disk0/zyy_dataset/Full/pokec.graph",
				          "/home/hnu/Disk0/zyy_dataset/pokec/part/p10",
			              "/home/hnu/Disk0/zyy_dataset/indochina/Index/index_",
			        1049, 123980, 10); //  Distri4hop 索引构建 的输入

	worker_finalize();
	

	return 0;
}
