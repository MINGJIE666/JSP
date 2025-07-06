#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <time.h>
#include <ctime>
#include <string>
#include <string.h>
#include <algorithm>
#include<math.h>
#include <iomanip>

#include<vector>

#define MAXN 52
#define MAXM 22

using namespace std;
string inFile = "ta01";//输入文件路径
int LB = 1231;
unsigned seed = 666;

string outFile = "output.txt";//输出文件路径
int endTime = 20;//程序运行时间
int tabuDepth = 200;//搜索深度
int searchDepth = 10e4;//搜索深度
clock_t startTime;//程序开始时间
double bestTime;

//下面是算法主要数据结构
int N, M; //工件数量N和机器数量M
int queue[MAXM * MAXN * 10], front, rear; //队列以及对应的头尾指针
int jm[MAXN][MAXM], jt[MAXN][MAXM];//仅用于读数据
int number[MAXM * MAXN][3];//1工件号，2工件的第几个工序
int job_time[MAXM * MAXN];//每个工序的加工时间
int job[MAXN * MAXM][5];//1工序开始时间，2工序结束时间，3所在的机器，4所在机器的第几个工序
int machine[MAXM][MAXN][3];//0机器M的第N个工序，1工序开始时间，2工序结束时间
int last_machine[MAXM];//每个机器最后完工的一个工序
int ans_makespan;//全局最优解
int this_makespan;//当前解
int HASH_LEN = 1E8;//哈希表的长度
int hash1[MAXM * MAXN][MAXM * MAXN], hash2[MAXM * MAXN][MAXM * MAXN], hash3[MAXM * MAXN][MAXM * MAXN];
int H1[100000000], H2[100000000], H3[100000000];
int hashObj1, hashObj2, hashObj3;
int r_change[MAXN], q_change[MAXN], critical[MAXM][MAXN];
int blocks[MAXM * MAXN * 2][3], q_save[MAXM * MAXN];
int min_makespan;
int choose_machine[MAXM * MAXN * MAXN], choose_u[MAXM * MAXN * MAXN], choose_v[MAXM * MAXN * MAXN], choose_forward[MAXM * MAXN * MAXN];
int choose[MAXM * MAXN * 100][4], ans_machine[MAXM][MAXN];
int start_line[MAXM], flag[MAXM][MAXN], flag_test[MAXM * MAXN];
int min_machine, min_u, min_v, min_forward;
int iter, bestIter;
int noImprove;//连续未更新的次数
int tabu1[MAXM * MAXN][MAXM * MAXN];
int tt = 10;//禁忌长度

int Qflag[MAXM * MAXN];
int start_flag[MAXM * MAXN];


void Read_Data();
void Init();
void Calculate_Q();
void Re_Coder();
void Jump_Pit(int jump_length, int per_flag);
bool If_Tabu(int u, int v, int fwd);
void Tabu_Search();
void show_solution();

inline int Random(int n) {//随机生成一个 1-n 的整数
	return rand() % n + 1;
}
inline int Mp(int i) {  //机器上的前序工序
	if (job[i][4] - 1 > 0)
		return machine[job[i][3]][job[i][4] - 1][0];
	else
		return 0;
}
inline int Ms(int i) {  //机器上的后序工序
	if (job[i][4] + 1 <= N)
		return machine[job[i][3]][job[i][4] + 1][0];
	else
		return 0;
}
inline int Jp(int i) {  //工件上的前序工序
	if (number[i][2] != 1)
		return (i - 1);
	else
		return 0;
}
inline int Js(int i) {   //工件上的后序工序
	if (number[i][2] != M)
		return (i + 1);
	else
		return 0;
}

void Read_Data() {//读数据
	int start_line[MAXN], temp_machine;
	memset(start_line, 0, sizeof(start_line));
	ifstream FIC;
	FIC.open(inFile);
	if (FIC.fail())
	{
		cout << "### Erreur open, File_Name ### " << inFile << endl;
		exit(0);
	}
	if (FIC.eof())
	{
		cout << "### Error open, File_Name ### " << inFile << endl;
		exit(0);
	}
	FIC >> N >> M;
	for (int loop1 = 1; loop1 <= N; loop1++) {
		for (int loop2 = 1; loop2 <= M; loop2++) {
			FIC >> jm[loop1][loop2] >> jt[loop1][loop2]; //记录工序 的执行机器和时间
			++jm[loop1][loop2];//机器从1开始
		}
	}
	FIC.close();
}

void Init() {//随机生成初始解
	int job_line[MAXN], machine_line[MAXM], temp_machine, random_job;
	memset(job_line, 0, sizeof(job_line));
	memset(machine_line, 0, sizeof(machine_line));
	memset(machine, 0, sizeof(machine));
	memset(job, 0, sizeof(job));
	memset(hash1, 0, sizeof(hash1));
	memset(hash2, 0, sizeof(hash2));
	memset(hash3, 0, sizeof(hash3));
	memset(H1, 0, sizeof(H1));
	memset(H2, 0, sizeof(H2));
	memset(H3, 0, sizeof(H3));

	for (int loop1 = 1; loop1 <= M * N; loop1++) {
		random_job = Random(N);
		while (job_line[random_job] >= M || random_job == 0) {
			random_job = (random_job + 1) % (N + 1);
		}
		temp_machine = jm[random_job][++job_line[random_job]];
		machine[temp_machine][++machine_line[temp_machine]][0] = (random_job - 1) * M + job_line[random_job];

		job[machine[temp_machine][machine_line[temp_machine]][0]][3] = temp_machine;
		job[machine[temp_machine][machine_line[temp_machine]][0]][4] = machine_line[temp_machine];
	}
	for (int loop1 = 1; loop1 <= N * M; loop1++) {  //number数组存job号和job内的序号 
		number[loop1][1] = ((loop1 - 1) / M) + 1;
		number[loop1][2] = ((loop1 - 1) % M) + 1;
	}
	job_time[0] = 0;
	for (int loop1 = 1; loop1 <= M * N; loop1++) {
		job_time[loop1] = jt[number[loop1][1]][number[loop1][2]];
	}

	//计算makespan
	int flag_ts[MAXM * MAXN];
	int temp_operation, jp, mp, js, ms;
	rear = front = 0;
	memset(flag_ts, 0, sizeof(flag_ts));
	//mj,如果一个工序是该机器的第一个工序则1，如果一个工序是一个工件的第一个工序则2，如果是3说明满足两个条件，先算这样的，加入queue
	for (int loop1 = 1; loop1 <= M; loop1++)
		flag_ts[machine[loop1][1][0]] = 1;
	for (int loop1 = 1; loop1 <= N; loop1++) {
		flag_ts[(loop1 - 1) * M + 1] += 2;
		if (flag_ts[(loop1 - 1) * M + 1] == 3)
			queue[rear++] = (loop1 - 1) * M + 1;
	}
	while (rear != front) {
		temp_operation = queue[front++];
		jp = Jp(temp_operation);
		mp = Mp(temp_operation);
		//if( job[0][2]!=0 )
		//	cout<<"wolegequ"<<endl;
		job[temp_operation][1] = machine[job[temp_operation][3]][job[temp_operation][4]][1] = max(job[mp][2], job[jp][2]);
		job[temp_operation][2] = machine[job[temp_operation][3]][job[temp_operation][4]][2] = job[temp_operation][1] + job_time[temp_operation];
		//flag[job[temp_operation][3]][job[temp_operation][4]]=-1;
		//start_line[job[temp_operation][3]]++;
		js = Js(temp_operation);
		ms = Ms(temp_operation);
		if (js != 0) {
			flag_ts[js] += 2;
			if (flag_ts[js] == 3)
				queue[rear++] = js;
		}
		if (ms != 0) {
			flag_ts[ms]++;
			if (flag_ts[ms] == 3)
				queue[rear++] = ms;
		}
	}

	int temp_max, loop;
	loop = 0;
	for (int loop1 = 1; loop1 <= M; loop1++) //计算makespan 
		if (loop1 == 1) {
			temp_max = 1;
			last_machine[++loop] = machine[loop1][N][0];
		}
		else if (machine[loop1][N][2] == machine[temp_max][N][2])
			last_machine[++loop] = machine[loop1][N][0];
		else if (machine[loop1][N][2] > machine[temp_max][N][2]) {
			temp_max = loop1;
			loop = 0;
			last_machine[++loop] = machine[loop1][N][0];
		}
	last_machine[0] = loop;//mj，last_machine用来记下相同结束时间的工序个数
	ans_makespan = this_makespan = machine[temp_max][N][2];//最后工序的总序号 

	//初始化哈希，给每条边赋随机值
	int temp_job1, temp_job2;
	memset(H1, 0, sizeof(H1));
	memset(H2, 0, sizeof(H2));
	memset(H3, 0, sizeof(H3));
	hashObj1 = hashObj2 = hashObj3 = 0;
	for (int loop1 = 1; loop1 <= M; loop1++) {
		for (int loop2 = 1; loop2 < N; loop2++) {
			temp_job1 = machine[loop1][loop2][0];
			for (int loop3 = loop2 + 1; loop3 <= N; loop3++) {
				temp_job2 = machine[loop1][loop3][0];
				hash1[temp_job1][temp_job2] = Random(1E5);
				hash1[temp_job2][temp_job1] = Random(1E5);
				hash2[temp_job1][temp_job2] = Random(1E5);
				hash2[temp_job2][temp_job1] = Random(1E5);
				hash3[temp_job1][temp_job2] = Random(1E5);
				hash3[temp_job2][temp_job1] = Random(1E5);
				if (loop3 == loop2 + 1) {//记录当前解的哈希值
					hashObj1 = (hashObj1 + hash1[temp_job1][temp_job2]) % HASH_LEN;
					hashObj2 = (hashObj2 + hash2[temp_job1][temp_job2]) % HASH_LEN;
					hashObj3 = (hashObj3 + hash3[temp_job1][temp_job2]) % HASH_LEN;
				}
			}
		}
	}

	H1[hashObj1] = 1;
	H2[hashObj2] = 1;
	H3[hashObj3] = 1;
}

void Critical_Path() {//寻找关键路径
	int operation, mp, jp, temp_machine, machine_number, loop3;
	front = rear = 0;
	memset(critical, 0, sizeof(critical));
	for (int loop1 = 1; loop1 <= last_machine[0]; loop1++)
		queue[rear++] = last_machine[loop1];
	while (front != rear) {                          //用队列实现广度优先，构造关键路径critical[][] 
		operation = queue[front++];
		mp = Mp(operation);
		jp = Jp(operation);
		temp_machine = job[operation][3];
		machine_number = job[operation][4];
		if (critical[temp_machine][machine_number] == 0) {
			critical[temp_machine][machine_number] = 1;
			if (mp != 0 && job[operation][1] == job[mp][2])
				queue[rear++] = mp;
			if (jp != 0 && job[operation][1] == job[jp][2])
				queue[rear++] = jp;
		}
	}

	loop3 = 0;
	for (int loop1 = 1; loop1 <= M; loop1++)       //构造关建块数组blocks[][] 
		for (int loop2 = 1; loop2 <= N; loop2++)
			if (critical[loop1][loop2] == 1 && loop2 <= N) {
				blocks[++loop3][0] = loop1;
				blocks[loop3][1] = loop2;//关键块为一个工序的情况 
				while (critical[loop1][loop2] == 1 && loop2 <= N) {
					loop2++;
				}
				blocks[loop3][2] = loop2 - 1;
			}
	blocks[0][0] = loop3;//存关键块的个数 
}

void Calculate_Q() {   //计算Q值
	int operation, rear, front;
	rear = front = q_save[0] = 0;
	memset(Qflag, 0, sizeof(flag));
	Qflag[0] = 1;
	//mj,记录每个机器最后一个工序，同时该工序也得是该工件的最后一个工序，存入queue
	for (int loop1 = 1; loop1 <= M; loop1++)
		Qflag[machine[loop1][N][0]]++;
	for (int loop1 = 1; loop1 <= N; loop1++) {
		Qflag[loop1 * M]++;
		if (Qflag[loop1 * M] == 2)
			queue[rear++] = loop1 * M;
	}
	while (front != rear) {
		operation = queue[front++];
		q_save[operation] = max(q_save[Js(operation)], q_save[Ms(operation)]) + job_time[operation];
		if (Mp(operation) != 0) {
			Qflag[Mp(operation)]++;
			if (Qflag[Mp(operation)] == 2)
				queue[rear++] = Mp(operation);
		}
		if (Jp(operation) != 0) {
			Qflag[Jp(operation)]++;
			if (Qflag[Jp(operation)] == 2)
				queue[rear++] = Jp(operation);
		}
	}

}

void Re_Coder() {  //调整数据结构
	int front, rear, temp_save;
	int f_machine_number = job[min_u][4];
	int e_machine_number = job[min_v][4];
	int temp_machine, temp_job1, temp_job2, temp_job3, temp_job4;
	int temp_h1, temp_h2, temp_h3;
	int u, v;
	u = min_u;
	v = min_v;
	temp_machine = job[u][3];
	temp_job1 = machine[temp_machine][job[u][4] - 1][0];
	temp_job2 = machine[temp_machine][job[u][4] + 1][0];
	temp_job3 = machine[temp_machine][job[v][4] - 1][0];
	temp_job4 = machine[temp_machine][job[v][4] + 1][0];
	if (min_forward == 1) {        //领域选择后  调整部分数据结构值,机器上的顺序
		temp_h1 = (hashObj1 + hash1[temp_job1][v] + hash1[v][u] + hash1[temp_job3][temp_job4] - hash1[temp_job1][u] - hash1[temp_job3][v] - hash1[v][temp_job4]) % HASH_LEN;
		temp_h2 = (hashObj2 + hash2[temp_job1][v] + hash2[v][u] + hash2[temp_job3][temp_job4] - hash2[temp_job1][u] - hash2[temp_job3][v] - hash2[v][temp_job4]) % HASH_LEN;
		temp_h3 = (hashObj3 + hash3[temp_job1][v] + hash3[v][u] + hash3[temp_job3][temp_job4] - hash3[temp_job1][u] - hash3[temp_job3][v] - hash3[v][temp_job4]) % HASH_LEN;
		temp_save = machine[min_machine][e_machine_number][0];
		for (int loop1 = e_machine_number; loop1 >= f_machine_number; loop1--) {
			if (loop1 == f_machine_number) {
				machine[min_machine][loop1][0] = temp_save;
				job[temp_save][4] = f_machine_number;
			}
			else {
				machine[min_machine][loop1][0] = machine[min_machine][loop1 - 1][0];
				job[machine[min_machine][loop1][0]][4] = loop1;
			}
		}
	}
	else {
		temp_h1 = (hashObj1 + hash1[temp_job1][temp_job2] + hash1[v][u] + hash1[u][temp_job4] - hash1[temp_job1][u] - hash1[u][temp_job2] - hash1[v][temp_job4]) % HASH_LEN;
		temp_h2 = (hashObj2 + hash2[temp_job1][temp_job2] + hash2[v][u] + hash2[u][temp_job4] - hash2[temp_job1][u] - hash2[u][temp_job2] - hash2[v][temp_job4]) % HASH_LEN;
		temp_h3 = (hashObj3 + hash3[temp_job1][temp_job2] + hash3[v][u] + hash3[u][temp_job4] - hash3[temp_job1][u] - hash3[u][temp_job2] - hash3[v][temp_job4]) % HASH_LEN;
		temp_save = machine[min_machine][f_machine_number][0];
		for (int loop1 = f_machine_number; loop1 <= e_machine_number; loop1++) {
			if (loop1 == e_machine_number) {
				machine[min_machine][loop1][0] = temp_save;
				job[temp_save][4] = e_machine_number;
			}
			else {

				machine[min_machine][loop1][0] = machine[min_machine][loop1 + 1][0];
				job[machine[min_machine][loop1][0]][4] = loop1;
			}
		}
	}
	H1[temp_h1] = H2[temp_h2] = H3[temp_h3] = 1;//评价过的解即标记
	hashObj1 = temp_h1;
	hashObj2 = temp_h2;
	hashObj3 = temp_h3;
	int  temp_operation;
	for (int loop1 = 0; loop1 <= M; loop1++)//初始化start_line
		start_line[loop1] = 100;
	start_line[min_machine] = f_machine_number;
	memset(start_flag, -1, sizeof(start_flag)); //算出start_line
	rear = front = 0;
	queue[rear++] = machine[min_machine][f_machine_number][0];
	start_flag[0] = 1;
	while (rear != front) {
		temp_operation = queue[front++];
		int js = Js(temp_operation);
		int ms = Ms(temp_operation);
		int temp_machine = job[temp_operation][3];
		start_line[temp_machine] = ((job[temp_operation][4] < start_line[temp_machine]) ? job[temp_operation][4] : start_line[temp_machine]);
		if (start_flag[js] != 1) {             //start_flag[]==1 ,则需要调整
			start_flag[js] = 1;
			queue[rear++] = js;
		}
		if (start_flag[ms] != 1) {
			start_flag[ms] = 1;
			queue[rear++] = ms;
		}
	}

	for (int loop1 = 1; loop1 <= M; loop1++)       //算出flag的值
		for (int loop2 = start_line[loop1]; loop2 <= N; loop2++) {
			flag[loop1][loop2] = 0;
			temp_operation = machine[loop1][loop2][0];
			int jp = Jp(temp_operation);
			int mp = Mp(temp_operation);
			if (job[mp][4] < start_line[job[mp][3]])//mp 不需要更改
				flag[loop1][loop2]++;
			if (job[jp][4] < start_line[job[jp][3]])//jp 不需要更改
				flag[loop1][loop2] += 2;
		}

	//clock_form4=clock();
	rear = front = 0;
	queue[rear++] = machine[min_machine][f_machine_number][0];
	while (rear != front) {
		temp_operation = queue[front++];
		int jp = Jp(temp_operation);
		int mp = Mp(temp_operation);

		job[temp_operation][1] = machine[job[temp_operation][3]][job[temp_operation][4]][1] = max(job[mp][2], job[jp][2]);
		job[temp_operation][2] = machine[job[temp_operation][3]][job[temp_operation][4]][2] = job[temp_operation][1] + job_time[temp_operation];
		int js = Js(temp_operation);
		int ms = Ms(temp_operation);
		if (js != 0) {
			flag[job[js][3]][job[js][4]] += 2;
			if (flag[job[js][3]][job[js][4]] == 3)
				queue[rear++] = js;
		}
		if (ms != 0) {
			flag[job[ms][3]][job[ms][4]]++;
			if (flag[job[ms][3]][job[ms][4]] == 3)
				queue[rear++] = ms;
		}

	}

	int temp_max, loop2;
	loop2 = 0;
	for (int loop1 = 1; loop1 <= M; loop1++) { //计算makespan 
		if (loop1 == 1) {
			temp_max = 1;
			last_machine[++loop2] = machine[1][N][0];
		}
		else if (machine[loop1][N][2] == machine[temp_max][N][2]) {
			last_machine[++loop2] = machine[loop1][N][0];
		}
		else if (machine[loop1][N][2] > machine[temp_max][N][2]) {
			temp_max = loop1;
			loop2 = 0;
			last_machine[++loop2] = machine[loop1][N][0];
		}
	}
	last_machine[0] = loop2;//记下相同结束时间的工序个数
	this_makespan = machine[temp_max][N][2];//最后工序的总序号
}

void Jump_Pit(int jump_length, int per_flag) {
	int u_machine_number, u_number, v_machine_number, v_number, job_number, temp_machine;
	int f_machine_number, e_machine_number, f_total_number, e_total_number;
	int jump_pit_cycle_number = jump_length;
	int loop;

	//实验发现不管什么情况都移动关键块会更好

	if (!per_flag) {//不是扰动的话，只动关键块
		while (jump_pit_cycle_number--) {
			Critical_Path();
			Calculate_Q();
			//memset(choose,0,sizeof(choose));
			loop = 0;
			for (int loop1 = 1; loop1 <= blocks[0][0]; loop1++) {
				e_machine_number = blocks[loop1][2];//结束节点的机器上的序号
				f_machine_number = blocks[loop1][1]; //开始节点的机器上的序号
				job_number = e_machine_number - f_machine_number + 1;//关键块中工序的个数
				if (job_number == 1) continue;
				temp_machine = blocks[loop1][0];
				f_total_number = machine[temp_machine][f_machine_number][0];//开始节点的总序号
				e_total_number = machine[temp_machine][e_machine_number][0];//结束节点的总序号
				for (int loop2 = 2; loop2 <= job_number; loop2++) {//一个关键块里面的工序序号
					u_machine_number = f_machine_number;//u的机器内序号
					u_number = machine[temp_machine][u_machine_number][0];//u的总序号
					v_machine_number = f_machine_number + loop2 - 1;//v的机器内序号
					v_number = machine[temp_machine][v_machine_number][0];//v的总序号
					if (q_save[v_number] >= q_save[Js(u_number)]) {//u放置v的可行解
						choose[++loop][0] = temp_machine;
						choose[loop][1] = u_number;
						choose[loop][2] = v_number;
						choose[loop][3] = 0;
					}
					if (job[u_number][2] >= job[Jp(v_number)][2]) {//v放置U 之前的可行解
						choose[++loop][0] = temp_machine;
						choose[loop][1] = u_number;
						choose[loop][2] = v_number;
						choose[loop][3] = 1;
					}
				}
				for (int loop2 = 1; loop2 <= job_number - 1; loop2++) {//blocks内序号
					u_machine_number = f_machine_number + loop2 - 1;//u的机器内序号
					u_number = machine[temp_machine][u_machine_number][0];//u的总序号
					v_machine_number = e_machine_number;
					v_number = machine[temp_machine][e_machine_number][0];//e_total_number
					if (q_save[v_number] >= q_save[Js(u_number)]) {//u放置V后面
						choose[++loop][0] = temp_machine;
						choose[loop][1] = u_number;
						choose[loop][2] = v_number;
						choose[loop][3] = 0;
					}
					if (job[u_number][2] >= job[Jp(v_number)][2]) {//v放置U前面
						choose[++loop][0] = temp_machine;
						choose[loop][1] = u_number;
						choose[loop][2] = v_number;
						choose[loop][3] = 1;
					}
				}
			}
			choose[0][0] = loop;//存可供选择的领域个数
			if (loop == 0)cout << "没有可选择的邻居解" << endl;
			int random_choose = Random(loop);//随机数 1----loop
			min_machine = choose[random_choose][0];
			min_u = choose[random_choose][1];
			min_v = choose[random_choose][2];
			min_forward = choose[random_choose][3];
			Re_Coder();
		}
	}
	else {//如果扰动就完全随机
		while (jump_pit_cycle_number--) {
			Critical_Path();
			Calculate_Q();
			loop = 0;
			for (int loop1 = 1; loop1 <= M; loop1++) {
				temp_machine = loop1;
				for (int loop2 = 1; loop2 < N; loop2++) {//一个关键块里面的工序序号 
					u_machine_number = loop2;//u的机器内序号 
					u_number = machine[temp_machine][u_machine_number][0];//u的总序号
					for (int loop3 = loop2 + 1; loop3 <= N; loop3++) {//一个关键块里面的工序序号				
						v_machine_number = loop3;//v的机器内序号 
						v_number = machine[temp_machine][v_machine_number][0];//v的总序号 
						if (q_save[v_number] >= q_save[Js(u_number)]) {//u放置v的可行解 
							choose[++loop][0] = temp_machine;
							choose[loop][1] = u_number;
							choose[loop][2] = v_number;
							choose[loop][3] = 0;
						}
						if (job[u_number][2] >= job[Jp(v_number)][2]) {//v放置U 之前的可行解 
							choose[++loop][0] = temp_machine;
							choose[loop][1] = u_number;
							choose[loop][2] = v_number;
							choose[loop][3] = 1;
						}
					}
				}
			}
			choose[0][0] = loop;//存可供选择的领域个数
			if (loop == 0)cout << "没有可选择的邻居解" << endl;
			int random_choose = Random(loop);//随机数 1----loop

			min_machine = choose[random_choose][0];
			min_u = choose[random_choose][1];
			min_v = choose[random_choose][2];
			min_forward = choose[random_choose][3];
			Re_Coder();
		}

	}

}

bool If_Tabu(int u, int v, int fwd) {//返回1表示解被禁忌，0表示不被禁忌，fwd0表示u放v后，fwd1表示v放u前

	int temp_machine, temp_job1, temp_job2, temp_job3, temp_job4;
	int temp_h1, temp_h2, temp_h3;
	temp_machine = job[u][3];
	temp_job1 = machine[temp_machine][job[u][4] - 1][0];
	temp_job2 = machine[temp_machine][job[u][4] + 1][0];
	temp_job3 = machine[temp_machine][job[v][4] - 1][0];
	temp_job4 = machine[temp_machine][job[v][4] + 1][0];
	if (!fwd) {//fwd0表示u放v后
		temp_h1 = (hashObj1 + hash1[temp_job1][temp_job2] + hash1[v][u] + hash1[u][temp_job4] - hash1[temp_job1][u] - hash1[u][temp_job2] - hash1[v][temp_job4]) % HASH_LEN;
		temp_h2 = (hashObj2 + hash2[temp_job1][temp_job2] + hash2[v][u] + hash2[u][temp_job4] - hash2[temp_job1][u] - hash2[u][temp_job2] - hash2[v][temp_job4]) % HASH_LEN;
		temp_h3 = (hashObj3 + hash3[temp_job1][temp_job2] + hash3[v][u] + hash3[u][temp_job4] - hash3[temp_job1][u] - hash3[u][temp_job2] - hash3[v][temp_job4]) % HASH_LEN;
	}
	else {//fwd1表示v放u前
		temp_h1 = (hashObj1 + hash1[temp_job1][v] + hash1[v][u] + hash1[temp_job3][temp_job4] - hash1[temp_job1][u] - hash1[temp_job3][v] - hash1[v][temp_job4]) % HASH_LEN;
		temp_h2 = (hashObj2 + hash2[temp_job1][v] + hash2[v][u] + hash2[temp_job3][temp_job4] - hash2[temp_job1][u] - hash2[temp_job3][v] - hash2[v][temp_job4]) % HASH_LEN;
		temp_h3 = (hashObj3 + hash3[temp_job1][v] + hash3[v][u] + hash3[temp_job3][temp_job4] - hash3[temp_job1][u] - hash3[temp_job3][v] - hash3[v][temp_job4]) % HASH_LEN;
	}

	if (H1[temp_h1] && H2[temp_h2] && H3[temp_h3]) {
			return true;
	}
	else {
		if (Random(100) <= 100) {//一定的概率禁忌
			H1[temp_h1] = H2[temp_h2] = H3[temp_h3] = 1;//评价过的解即标记 
		}
		return false;//不被禁忌
	}
}

void Tabu_Search() {//禁忌搜索
	int temp_makespan, temp_machine, u_machine_number, u_number, v_machine_number, v_number, w_number, job_number;
	int choose_number, f_machine_number, e_machine_number;
	int end_number;//关键块开头到机器结尾的数量
	int start_number;//关键块结尾到机器开始的数量

	min_makespan = HASH_LEN, choose_number = 0;
	Critical_Path();//计算关键路径

	Calculate_Q();//计算Q值;

	//N6邻域
	for (int loop1 = 1; loop1 <= blocks[0][0]; loop1++) {
		e_machine_number = blocks[loop1][2];//结束节点的机器上的序号
		f_machine_number = blocks[loop1][1]; //开始节点的机器上的序号 
		job_number = e_machine_number - f_machine_number + 1;//关键块中工序的个数 
		if (job_number == 1) continue;
		temp_machine = blocks[loop1][0];
		for (int loop2 = 2; loop2 <= job_number; loop2++) {//一个关键块里面的V的相对位置 
			u_machine_number = f_machine_number;//u的机器内序号 
			u_number = machine[temp_machine][u_machine_number][0];//u的总序号
			v_machine_number = f_machine_number + loop2 - 1;//v的机器内序号 
			v_number = machine[temp_machine][v_machine_number][0];//v的总序号 
			if (q_save[v_number] >= q_save[Js(u_number)] && Jp(u_number) != 0 && !If_Tabu(u_number, v_number, 0)) {//u放置v的可行解 
				for (int loop3 = 2; loop3 <= loop2 + 1; loop3++) {//(0,w)
					w_number = machine[temp_machine][f_machine_number + loop3 - 1][0];//w的总序号 
					if (loop3 == 2 && Mp(u_number) == 0)//如果u是机器上第一个工序的话  L1  1
						r_change[loop3] = job[Jp(w_number)][2];
					else if (loop3 == 2)//L1  2
						r_change[loop3] = max(job[Jp(w_number)][2], job[Mp(u_number)][2]);
					else if (loop3 == loop2 + 1)//u
						r_change[loop3] = max(job[Jp(u_number)][2], r_change[loop3 - 1] + job_time[v_number]);//新图中的(0,u)
					else //w
						r_change[loop3] = max(job[Jp(w_number)][2], r_change[loop3 - 1] + job_time[Mp(w_number)]);
				}
				for (int loop3 = loop2 + 1; loop3 >= 2; loop3--) {//(w,n)
					w_number = machine[temp_machine][f_machine_number + loop3 - 1][0];//w的总序号
					if (loop3 == loop2 + 1 && Ms(v_number) == 0)//如果v是机器上的最后一个工序的话   U 1
						q_change[loop3] = job_time[u_number] + q_save[Js(u_number)];
					else if (loop3 == loop2 + 1)//U 2
						q_change[loop3] = job_time[u_number] + max(q_save[Js(u_number)], q_save[Ms(v_number)]);
					else if (loop3 == loop2)//v
						q_change[loop3] = job_time[v_number] + max(q_save[Js(v_number)], q_change[loop3 + 1]);
					else //w
						q_change[loop3] = job_time[w_number] + max(q_save[Js(w_number)], q_change[loop3 + 1]);
				}

				temp_makespan = 0;
				for (int loop3 = 2; loop3 <= loop2 + 1; loop3++)//整合(0,n)的最大值 即估计的makespan 
					if (r_change[loop3] + q_change[loop3] > temp_makespan)
						temp_makespan = r_change[loop3] + q_change[loop3];

				if (temp_makespan <= min_makespan) {
					if (temp_makespan < min_makespan) {
						min_makespan = temp_makespan;
						choose_number = 1;
						choose_machine[choose_number] = temp_machine;
						choose_u[choose_number] = u_number;
						choose_v[choose_number] = v_number;
						choose_forward[choose_number] = 0;
						//choose_forward[0]=0;
					}
					else if (temp_makespan == min_makespan) {
						choose_number++;
						choose_machine[choose_number] = temp_machine;
						choose_u[choose_number] = u_number;
						choose_v[choose_number] = v_number;
						choose_forward[choose_number] = 0;
						//choose_forward[0]=0;
					}
				}

			}

			if (job[u_number][2] >= job[Jp(v_number)][2] && Jp(u_number) != 0 && !If_Tabu(u_number, v_number, 1)) {//v放置U之前的可行解 
				for (int loop3 = 0; loop3 <= loop2 - 1; loop3++) {//(0,w)
					w_number = machine[temp_machine][f_machine_number + loop3 - 1][0];
					if (loop3 == 0 && Mp(u_number) == 0)//v 1     f_machine_number-1==0
						r_change[loop3] = job[Jp(v_number)][2];
					else if (loop3 == 0)//v 2
						r_change[loop3] = max(job[Jp(v_number)][2], job[Mp(u_number)][2]);
					else if (loop3 == 1)//u
						r_change[loop3] = max(job[Jp(u_number)][2], r_change[loop3 - 1] + job_time[v_number]);
					else //w
						r_change[loop3] = max(job[Jp(w_number)][2], r_change[loop3 - 1] + job_time[Mp(w_number)]);
				}

				for (int loop3 = loop2 - 1; loop3 >= 0; loop3--) {//(w,n)
					w_number = machine[temp_machine][f_machine_number + loop3 - 1][0];
					if (loop3 == loop2 - 1 && Ms(v_number) == 0)//Lk 1
						q_change[loop3] = job_time[w_number] + q_save[Js(w_number)];
					else if (loop3 == loop2 - 1)//Lk 2
						q_change[loop3] = job_time[w_number] + max(q_save[Js(w_number)], q_save[Ms(v_number)]);
					else if (loop3 == 0)//v
						q_change[loop3] = job_time[v_number] + max(q_save[Js(v_number)], q_change[loop3 + 1]);
					else //w
						q_change[loop3] = job_time[w_number] + max(q_save[Js(w_number)], q_change[loop3 + 1]);
				}

				temp_makespan = 0;
				for (int loop3 = 0; loop3 <= loop2 - 1; loop3++)//整合(0,n)的最大值 即估计的makespan 
					if (r_change[loop3] + q_change[loop3] > temp_makespan)
						temp_makespan = r_change[loop3] + q_change[loop3];

				if (temp_makespan <= min_makespan) {
					if (temp_makespan < min_makespan) {
						min_makespan = temp_makespan;
						choose_number = 1;
						choose_machine[choose_number] = temp_machine;
						choose_u[choose_number] = u_number;
						choose_v[choose_number] = v_number;
						choose_forward[choose_number] = 1;
						//choose_forward[0]=0;
					}
					else if (temp_makespan == min_makespan) {
						choose_number++;
						choose_machine[choose_number] = temp_machine;
						choose_u[choose_number] = u_number;
						choose_v[choose_number] = v_number;
						choose_forward[choose_number] = 1;
						//choose_forward[0]=0;
					}
				}
			}
		}

		for (int loop2 = 1; loop2 <= job_number - 1; loop2++) {//blocks内序号 
			u_machine_number = f_machine_number + loop2 - 1;//u的机器内序号 
			u_number = machine[temp_machine][u_machine_number][0];//u的总序号
			v_machine_number = e_machine_number;
			v_number = machine[temp_machine][e_machine_number][0];//e_total_number
			if (q_save[v_number] >= q_save[Js(u_number)] && Js(v_number) != 0 && !If_Tabu(u_number, v_number, 0)) {//u放置v的可行解 
				for (int loop3 = loop2 + 1; loop3 <= job_number + 1; loop3++) {//(0,w)    e_blocks_number=job_number
					w_number = machine[temp_machine][f_machine_number + loop3 - 1][0];//w的总序号 
					if (loop3 == loop2 + 1 && Mp(u_number) == 0)//L1  1
						r_change[loop3] = job[Jp(w_number)][2];
					else if (loop3 == loop2 + 1)//L1  2
						r_change[loop3] = max(job[Jp(w_number)][2], job[Mp(u_number)][2]);
					else if (loop3 == job_number + 1)//u
						r_change[loop3] = max(job[Jp(u_number)][2], r_change[loop3 - 1] + job_time[v_number]);
					else //w
						r_change[loop3] = max(job[Jp(w_number)][2], r_change[loop3 - 1] + job_time[Mp(w_number)]);
				}
				for (int loop3 = job_number + 1; loop3 >= loop2 + 1; loop3--) {// (w,n)
					int w_number = machine[temp_machine][f_machine_number + loop3 - 1][0];//w的总序号 
					if (loop3 == job_number + 1 && Ms(v_number) == 0)//u  1
						q_change[loop3] = job_time[u_number] + q_save[Js(u_number)];
					else if (loop3 == job_number + 1)//u  2
						q_change[loop3] = job_time[u_number] + max(q_save[Js(u_number)], q_save[Ms(v_number)]);
					else if (loop3 == job_number)//v
						q_change[loop3] = job_time[v_number] + max(q_save[Js(v_number)], q_change[loop3 + 1]);
					else //w
						q_change[loop3] = job_time[w_number] + max(q_save[Js(w_number)], q_change[loop3 + 1]);
				}
				temp_makespan = 0;
				for (int loop3 = loop2 + 1; loop3 <= job_number + 1; loop3++)//整合(0,n)的最大值 即估计的makespan 
					if (r_change[loop3] + q_change[loop3] > temp_makespan)
						temp_makespan = r_change[loop3] + q_change[loop3];

				if (temp_makespan <= min_makespan) {
					if (temp_makespan < min_makespan) {
						min_makespan = temp_makespan;
						choose_number = 1;
						choose_machine[choose_number] = temp_machine;
						choose_u[choose_number] = u_number;
						choose_v[choose_number] = v_number;
						choose_forward[choose_number] = 0;
						//	choose_forward[0]=0;
					}
					else if (temp_makespan == min_makespan) {
						choose_number++;
						choose_machine[choose_number] = temp_machine;
						choose_u[choose_number] = u_number;
						choose_v[choose_number] = v_number;
						choose_forward[choose_number] = 0;
						//	choose_forward[0]=0;
					}
				}
			
			}
			if (job[u_number][2] >= job[Jp(v_number)][2] && Js(v_number) != 0 && !If_Tabu(u_number, v_number, 1)) {//v放置U前面 
				for (int loop3 = loop2 - 1; loop3 <= job_number - 1; loop3++) {//(0,w)   
					w_number = machine[temp_machine][f_machine_number + loop3 - 1][0];//w的总序号    
					if (loop3 == (loop2 - 1) && Mp(u_number) == 0)   //v  1
						r_change[loop3] = job[Jp(v_number)][2];
					else if (loop3 == loop2 - 1)//v 2
						r_change[loop3] = max(job[Jp(v_number)][2], job[Mp(u_number)][2]);
					else if (loop3 == loop2)//u
						r_change[loop3] = max(job[Jp(u_number)][2], r_change[loop3 - 1] + job_time[v_number]);
					else //w
						r_change[loop3] = max(job[Jp(w_number)][2], r_change[loop3 - 1] + job_time[Mp(w_number)]);
				}
				for (int loop3 = job_number - 1; loop3 >= loop2 - 1; loop3--) {//(w,n)
					w_number = machine[temp_machine][f_machine_number + loop3 - 1][0];//w的总序号
					if (loop3 == (job_number - 1) && Ms(v_number) == 0)//lk 1
						q_change[loop3] = job_time[w_number] + q_save[Js(w_number)];
					else if (loop3 == job_number - 1)//lk 2
						q_change[loop3] = job_time[w_number] + max(q_save[Js(w_number)], q_save[Ms(v_number)]);
					else if (loop3 == loop2 - 1)//v
						q_change[loop3] = job_time[v_number] + max(q_save[Js(v_number)], q_change[loop3 + 1]);
					else //w
						q_change[loop3] = job_time[w_number] + max(q_save[Js(w_number)], q_change[loop3 + 1]);
				}
				temp_makespan = 0;
				for (int loop3 = loop2 - 1; loop3 <= job_number - 1; loop3++)//整合(0,n)的最大值 即估计的makespan 
					if (r_change[loop3] + q_change[loop3] > temp_makespan)
						temp_makespan = r_change[loop3] + q_change[loop3];

				////
				if (temp_makespan <= min_makespan) {
					if (temp_makespan < min_makespan) {
						min_makespan = temp_makespan;
						choose_number = 1;
						choose_machine[choose_number] = temp_machine;
						choose_u[choose_number] = u_number;
						choose_v[choose_number] = v_number;
						choose_forward[choose_number] = 1;
						//	choose_forward[0]=0;
					}
					else if (temp_makespan == min_makespan) {
						choose_number++;
						choose_machine[choose_number] = temp_machine;
						choose_u[choose_number] = u_number;
						choose_v[choose_number] = v_number;
						choose_forward[choose_number] = 1;
						//choose_forward[0]=0;
					}
				}
			
			}
		}
	}

	//将u移动到关键块后面
#if 1 
	for (int loop1 = 1; loop1 <= blocks[0][0]; loop1++) {
		e_machine_number = N;//结束节点的机器上的序号
		f_machine_number = blocks[loop1][1]; //开始节点的机器上的序号 
		job_number = e_machine_number - f_machine_number + 1;//关键块中工序的个数 
		if (job_number == 1) continue;
		temp_machine = blocks[loop1][0];
		for (int loop2 = blocks[loop1][2]+1; loop2 <= job_number; loop2++) {//一个关键块里面的V的相对位置 
			u_machine_number = f_machine_number;//u的机器内序号 
			u_number = machine[temp_machine][u_machine_number][0];//u的总序号
			v_machine_number = f_machine_number + loop2 - 1;//v的机器内序号 
			v_number = machine[temp_machine][v_machine_number][0];//v的总序号 
			if (q_save[v_number] >= q_save[Js(u_number)] && Jp(u_number) != 0 && !If_Tabu(u_number, v_number, 0)) {//u放置v的可行解 
				for (int loop3 = 2; loop3 <= loop2 + 1; loop3++) {//(0,w)
					w_number = machine[temp_machine][f_machine_number + loop3 - 1][0];//w的总序号 
					if (loop3 == 2 && Mp(u_number) == 0)//如果u是机器上第一个工序的话  L1  1
						r_change[loop3] = job[Jp(w_number)][2];
					else if (loop3 == 2)//L1  2
						r_change[loop3] = max(job[Jp(w_number)][2], job[Mp(u_number)][2]);
					else if (loop3 == loop2 + 1)//u
						r_change[loop3] = max(job[Jp(u_number)][2], r_change[loop3 - 1] + job_time[v_number]);//新图中的(0,u)
					else //w
						r_change[loop3] = max(job[Jp(w_number)][2], r_change[loop3 - 1] + job_time[Mp(w_number)]);
				}
				for (int loop3 = loop2 + 1; loop3 >= 2; loop3--) {//(w,n)
					w_number = machine[temp_machine][f_machine_number + loop3 - 1][0];//w的总序号
					if (loop3 == loop2 + 1 && Ms(v_number) == 0)//如果v是机器上的最后一个工序的话   U 1
						q_change[loop3] = job_time[u_number] + q_save[Js(u_number)];
					else if (loop3 == loop2 + 1)//U 2
						q_change[loop3] = job_time[u_number] + max(q_save[Js(u_number)], q_save[Ms(v_number)]);
					else if (loop3 == loop2)//v
						q_change[loop3] = job_time[v_number] + max(q_save[Js(v_number)], q_change[loop3 + 1]);
					else //w
						q_change[loop3] = job_time[w_number] + max(q_save[Js(w_number)], q_change[loop3 + 1]);
				}

				temp_makespan = 0;
				for (int loop3 = 2; loop3 <= loop2 + 1; loop3++)//整合(0,n)的最大值 即估计的makespan 
					if (r_change[loop3] + q_change[loop3] > temp_makespan)
						temp_makespan = r_change[loop3] + q_change[loop3];

				if (temp_makespan <= min_makespan) {
					if (temp_makespan < min_makespan) {
						min_makespan = temp_makespan;
						choose_number = 1;
						choose_machine[choose_number] = temp_machine;
						choose_u[choose_number] = u_number;
						choose_v[choose_number] = v_number;
						choose_forward[choose_number] = 0;
						//choose_forward[0]=0;
						//cout << "N3" << endl;

					}
					else if (temp_makespan == min_makespan) {
						choose_number++;
						choose_machine[choose_number] = temp_machine;
						choose_u[choose_number] = u_number;
						choose_v[choose_number] = v_number;
						choose_forward[choose_number] = 0;
						//choose_forward[0]=0;
					}
				}

			}
		}
	}
#endif

	//将v移动到关键块前面
#if 1 
	for (int loop1 = 1; loop1 <= blocks[0][0]; loop1++) {
		e_machine_number = blocks[loop1][2];//结束节点的机器上的序号
		f_machine_number = 1; //开始节点的机器上的序号 
		job_number = e_machine_number - f_machine_number + 1;//关键块中工序的个数 
		if (job_number == 1) continue;
		temp_machine = blocks[loop1][0];

		for (int loop2 = 1; loop2 <= job_number - 1; loop2++) {//blocks内序号 
			u_machine_number = f_machine_number + loop2 - 1;//u的机器内序号 
			if (u_machine_number >= blocks[loop1][1])break;
			u_number = machine[temp_machine][u_machine_number][0];//u的总序号
			v_machine_number = e_machine_number;
			v_number = machine[temp_machine][e_machine_number][0];//e_total_number
			if (job[u_number][2] >= job[Jp(v_number)][2] && Js(v_number) != 0 && !If_Tabu(u_number, v_number, 1)) {//v放置U前面 
				for (int loop3 = loop2 - 1; loop3 <= job_number - 1; loop3++) {//(0,w)   
					w_number = machine[temp_machine][f_machine_number + loop3 - 1][0];//w的总序号    
					if (loop3 == (loop2 - 1) && Mp(u_number) == 0)   //v  1
						r_change[loop3] = job[Jp(v_number)][2];
					else if (loop3 == loop2 - 1)//v 2
						r_change[loop3] = max(job[Jp(v_number)][2], job[Mp(u_number)][2]);
					else if (loop3 == loop2)//u
						r_change[loop3] = max(job[Jp(u_number)][2], r_change[loop3 - 1] + job_time[v_number]);
					else //w
						r_change[loop3] = max(job[Jp(w_number)][2], r_change[loop3 - 1] + job_time[Mp(w_number)]);
				}
				for (int loop3 = job_number - 1; loop3 >= loop2 - 1; loop3--) {//(w,n)
					w_number = machine[temp_machine][f_machine_number + loop3 - 1][0];//w的总序号
					if (loop3 == (job_number - 1) && Ms(v_number) == 0)//lk 1
						q_change[loop3] = job_time[w_number] + q_save[Js(w_number)];
					else if (loop3 == job_number - 1)//lk 2
						q_change[loop3] = job_time[w_number] + max(q_save[Js(w_number)], q_save[Ms(v_number)]);
					else if (loop3 == loop2 - 1)//v
						q_change[loop3] = job_time[v_number] + max(q_save[Js(v_number)], q_change[loop3 + 1]);
					else //w
						q_change[loop3] = job_time[w_number] + max(q_save[Js(w_number)], q_change[loop3 + 1]);
				}
				temp_makespan = 0;
				for (int loop3 = loop2 - 1; loop3 <= job_number - 1; loop3++)//整合(0,n)的最大值 即估计的makespan 
					if (r_change[loop3] + q_change[loop3] > temp_makespan)
						temp_makespan = r_change[loop3] + q_change[loop3];

				////
				if (temp_makespan <= min_makespan) {
					if (temp_makespan < min_makespan) {
						min_makespan = temp_makespan;
						choose_number = 1;
						choose_machine[choose_number] = temp_machine;
						choose_u[choose_number] = u_number;
						choose_v[choose_number] = v_number;
						choose_forward[choose_number] = 1;
						//	choose_forward[0]=0;
						//cout << "N4" << endl;
					}
					else if (temp_makespan == min_makespan) {
						choose_number++;
						choose_machine[choose_number] = temp_machine;
						choose_u[choose_number] = u_number;
						choose_v[choose_number] = v_number;
						choose_forward[choose_number] = 1;
						//choose_forward[0]=0;
					}
				}


			}
		}




	}
#endif

	if (choose_number == 0) {
		Jump_Pit(2, 0);
	}
	else {
		choose_number = Random(choose_number);
		min_machine = choose_machine[choose_number];
		min_u = choose_u[choose_number];
		min_v = choose_v[choose_number];
		min_forward = choose_forward[choose_number];
		Re_Coder();
	}


}

void show_solution() {//展示当前解

	ofstream fout(outFile);
	//fout << inFile << endl << endl;
	fout << "************gantt chart:************" << endl;
	for (int loop1 = 1; loop1 <= M; loop1++) {
		for (int loop2 = 1; loop2 <= N; loop2++) {
			fout << setw(2) << number[machine[loop1][loop2][0]][1] << "," << setw(2) << number[machine[loop1][loop2][0]][2]
				<< ": [" << left << setw(5) << machine[loop1][loop2][1] << " "
				<< right << setw(5) << machine[loop1][loop2][2] << "]	";
		}
		fout << endl;
	}
	fout << endl;


	for (int loop1 = 1; loop1 <= M; loop1++) {
		for (int loop2 = 1; loop2 <= N; loop2++) {
			cout << setw(5) << machine[loop1][loop2][0] << "	";
		}
		cout << endl;
	}
	cout << "hashObj1=" << hashObj1 << endl;
}

int main(int argc, char** argv)
{
	/*if (argc == 4){
		inFile = argv[1];
		seed = (unsigned)atoi(argv[2]);
		LB = (int)atoi(argv[3]);
	}
	else{
		cout << "Error : the user should input infile to run the program." << endl;
		exit(0);
	}*/
	
	startTime = clock();
	iter = 0;
	//seed = (unsigned)time(NULL);
	srand(seed);
	Read_Data();
	Init();
	cout << fixed << setprecision(3);


	while ((clock()-startTime)/CLOCKS_PER_SEC < endTime && ans_makespan>LB) {//五千万次迭代（高亮N8）
	//while (time(NULL)-t0<endTime) {
		Tabu_Search();
		if (this_makespan < ans_makespan) {//更新最优解
			bestTime = (double)(clock() - startTime) / CLOCKS_PER_SEC;
			noImprove = 0;
			ans_makespan = this_makespan;
			bestIter = iter;
			cout << "Find the better solution: " << left << setw(6) << ans_makespan << bestTime << "s iter=" << iter << endl;
		}
		else {
			noImprove++;
			if (noImprove % tabuDepth == 0) {//此处要进行扰动		
				Jump_Pit(2, 0);
				if (noImprove == searchDepth) {//如果连续5次扰动，重置哈希
					Jump_Pit(10, 0);
					noImprove = 0;
					memset(H1, 0, sizeof(H1));
					memset(H2, 0, sizeof(H2));
					memset(H3, 0, sizeof(H3));
				}
			}
		}

		//if(iter%10000==0)cout << "The best found solution: " << left << setw(6) << ans_makespan << bestTime << "s iter=" << bestIter << endl;

		++iter;
	}
	//cout << "ALL FINISH" << endl;
	cout << inFile << "	" << ans_makespan << "	" << bestTime << "	"<< bestIter<<"	" << seed << endl;
	return 0;
}

