#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<vector>
#include<time.h>
#define N_limit 80000000
using namespace std;
const int multiplier=500;
const int row=2000;
const int col=2000;
const int G_row=200;
const int G_col=200;
const double pressure_h2=0.1;
const double saturation=0.6358;
const double pressure_ch4=10;
const double kT=0.112;  // 1300 K
int num_evt[38]={0};
int num_detach_h[2]={0};
int num_attach_h[2]={0};
double num_diff[7]={0};
int mesh[row][col][2]={0};
int attach_num[8]={0};
int attach_num_h[2]={0};
int species_num[8]={0};

//add 09-12
double  S=row*col-(G_row+2)*(G_col+2);  //area of bare cu surface
double N_H=0;  //number of H adatoms
double D_H=0;  //density of H adatoms
int num_onsite[9]={0};  //0-5: number of species on graphene edges. C,CH,CH2,C2,C2H,C2H2
int num_det[9]={0};    //number of detach
double prob_onsite[9]={1.50,1.57,1.08,2.08,2.29,2.31,3.55,2.90,3.30}; //barriers of detachment
double prob_each_h[2]={1.20,1.59};//barriers of CH and C2 attachment onto h-saturated 

double prob_each[9][14]={
//	  0    1    2     3    4    5   6    7    8    9    10   11   12  13
//   diff  H    C    CH   CH2  CH3  C2  C2H  C2H2  Gr   e1   e2   e3  e4
	{0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00},
	{0.13,0.80,0.79,0.65,0.68,0.69,0.72,0.85,0.00,0.58,0.00,0.00,0.00,0.80},//H 
	{0.50,0.79,0.25,1.27,0.00,0.00,0.00,0.00,0.00,1.27,0.00,0.00,0.00,0.79},//C 
	{0.10,0.65,1.27,0.14,0.00,0.00,0.00,0.00,0.00,0.44,1.32,0.00,0.00,0.65},//CH    e1:CH->C+H
	{0.18,0.68,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.19,1.21,0.00,0.00,0.68},//CH2   e1:CH2->CH+H
	{0.00,0.69,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,1.39,0.00,0.00,0.69},//CH3   e1:CH3->CH2+H
	{0.49,0.72,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.58,2.65,0.00,0.00,0.72},//C2    e1:C2-> C+C
	{0.32,0.85,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.56,1.38,2.95,0.00,0.85},//C2H   e1:C2H->C2+H; e2:C2H->C+CH
	{0.44,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,1.05,1.64,2.28,1.83,0.00} //C2H2  e1:C2H2->C2H+H; e2:C2H2->CH+CH; e3: C2H2 desorp
};
int sp_merge[8][8]={\
//	1	2	3	4	5	6	7	8
//	H	C	CH	CH2	CH3	C2	C2H	C2H2
	{-1,	3,	4,	5,	-1,	7,	8,	0},//H
	{3,	6,	7,	0,	0,	0,	0,	0},//C
	{4,	7,	8,	0,	0,	0,	0,	0},//CH
	{5,	0,	0,	0,	0,	0,	0,	0},//CH2
	{-1,	0,	0,	0,	0,	0,	0,	0},//CH3
	{7,	0,	0,	0,	0,	0,	0,	0},//C2
	{8,	0,	0,	0,	0,	0,	0,	0},//C2H
	{0,	0,	0,	0,	0,	0,	0,	0}//C2H2
};

//revised 0913
//-1 do nothing
//-2 add H adatom, actually only N_H++
//-4 N_H--
int sp_e[8][4][2]={\
	{{0,0},{0,0},{0,0},{-4,-4}},  //H
	{{0,0},{0,0},{0,0},{3,-4}},  //C
	{{-2,2},{0,0},{0,0},{4,-4}},  //CH
	{{-2,3},{0,0},{0,0},{5,-4}},  //CH2
	{{-2,4},{0,0},{0,0},{-4,-1}},  //CH3
	{{2,2},{0,0},{0,0},{7,-4}},  //C2
	{{-2,6},{2,3},{0,0},{8,-4}},  //C2H
	{{-2,7},{3,3},{-1,-1},{0,0}}  //C2H2
};

int evt_emerge[8][8]={\
//	1	2	3	4	5	6	7	8
//	H	C	CH	CH2	CH3	C2	C2H	C2H2
	{2,	16,	15,	14,	13,	10,	11,	-1},	//H
	{16,	7,	8,	-1,	-1,	-1,	-1,	-1},	//C
	{15,	8,	9,	-1,	-1,	-1,	-1,	-1},	//CH
	{14,	-1,	-1,	-1,	-1,	-1,	-1,	-1},	//CH2
	{13,	-1,	-1,	-1,	-1,	-1,	-1,	-1},	//CH3
	{10,	-1,	-1,	-1,	-1,	-1,	-1,	-1},	//C2
	{11,	-1,	-1,	-1,	-1,	-1,	-1,	-1},	//C2H
	{-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1}	//C2H2
};

int evt_spec[8][4]={\
	{-1,	-1,	-1,	2},	//H
	{-1,	-1,	-1,	16},	//C
	{6,	-1,	-1,	15},	//CH
	{5,	-1,	-1,	14},	//CH2
	{4,	-1,	-1,	13},	//CH3
	{17,	-1,	-1,	10},	//C2
	{20,	18,	-1,	11},	//C2H
	{21,	19,	30,	-1}	//C2H2
};

int p[3][2];
int q[3][2];

void err(char* erro);
void init();
int get_nb(int x,int y,int xy,int i);
void find_evt();
void adsorption(int sp1,int sp2);
void diffusion(int sp,int s_n,int evt);
void merge(int sp1,int s_n,int evt);
void spec_del(int sp,int s_n);
void creat(int sp,int x,int y);
void spec_e(int sp,int s_n,int evt);
void update_evt_nb(int x,int y);
void statistic(int n,int N);
int get_sp(int s_x,int s_y);
void update_H();
void detach(int sp);

class species{
	public:
	species(int sp,int x,int y){
		spec=sp;
		s_x=x;
		s_y=y;
		prob[3]=prob_each[spec][10]; //desorp or decomp events 
		prob[4]=prob_each[spec][11];
		prob[5]=prob_each[spec][12];
		prob[6]=D_H*prob_each[spec][13]; //or coalescent events
		check_event();
	}
	void check_event(){
		double r9;
		for(int i=0;i<3;i++){
			nb[i][0]=get_nb(s_x,s_y,0,i);
			nb[i][1]=get_nb(s_x,s_y,1,i);
			event[i]=get_sp(nb[i][0],nb[i][1]);
			if(event[i]==100) prob[i]=0;
			//2017 added for h
			else if(event[i]==9) {
				//printf("event-9 spec: %d\n",spec);
				r9=(double)rand()/RAND_MAX;
				if(r9>saturation) prob[i]=prob_each[spec][event[i]];
				else {
					event[i]=10;
					if(spec==3)  prob[i]=prob_each_h[0];
					else if(spec==6)  prob[i]=prob_each_h[1];
					else prob[i]=0;
					//printf("spec  %d\t%f\n",spec,prob[i]);
				}
			}
			else prob[i]=prob_each[spec][event[i]];
		}
		prob_s=0;
		for(int i=0;i<7;i++) prob_s+=prob[i];
	}
	int update(){
		int i;
		double evt_sum=0;
		double evt_add=0;
		double evt_rand;
		for(i=0;i<7;i++)
			if(prob[i]!=0) evt_sum+=prob[i];
		evt_rand=evt_sum*(double)rand()/RAND_MAX;
		for(i=0;i<7;i++){
			if(prob[i]!=0) evt_add+=prob[i];
			if(evt_add>=evt_rand){
				return i;
			}
		}
		err("spec.update()");
		return 0;
	}
	void upd_H(){
		prob_s-=prob[6];
		if(spec!=4) prob[6]=D_H*prob_each[spec][13]; //or coalescent events
		else prob[6]=D_H*prob_each[spec][13]; //or coalescent events
		prob_s+=prob[6];
	}

	double prob[7],prob_s;
	int event[3],nb[3][2];
	int spec,s_x,s_y;
};	

vector<vector<species> > spec_list(8);

int main(){
	init();
	int init_h=0,h;
	for(int i=1;i<8;i++) init_h+=species_num[i];
	for(int m=0;m<multiplier;m++){
		for(int i=0;i<N_limit;i++){
			h=init_h;
			find_evt();
			init_h=0;
			for(int j=1;j<8;j++) init_h+=species_num[j];
			if(init_h!=h) statistic(i,m);
			if(i%40000000==0) statistic(-1,0);
		}
	}
	statistic(-1,0);
	for(int i=0;i<38;i++)
		printf("num_evt	%d:\t%d\n",i,num_evt[i]);
	for(int i=0;i<2;i++)
		printf("num_attach-h	%d:\t%d\n",i,num_attach_h[i]);
	for(int i=0;i<2;i++)
		printf("num_detach-h	%d:\t%d\n",i,num_detach_h[i]);
	for(int i=0;i<7;i++)
		printf("num_diff	%d:\t%f\n",i,num_diff[i]);
	return 0;
}

void err(char* erro){
	printf("%s error!!!\n",erro);
	scanf("%s",&erro);
}

void init(){
	srand((unsigned)time(NULL));
	for(int i=1;i<9;i++)
		for(int j=0;j<14;j++)
			if(prob_each[i][j]!=0)
				prob_each[i][j]=2.71*pow(10,13)*exp(-prob_each[i][j]/0.112);
	for(int i=1;i<9;i++)
		for(int j=1;j<9;j++)
			if(prob_each[i][j]!=0)
				prob_each[i][j]=prob_each[i][j]/2;
	for(int i=0;i<9;i++){ prob_onsite[i]=2.71*pow(10,13)*exp(-prob_onsite[i]/0.112);
		printf("prob_onsite: %d\t%f\n",i,prob_onsite[i]);}
	for(int i=0;i<2;i++){ prob_each_h[i]=2.71*pow(10,13)*exp(-prob_each_h[i]/0.112);
		printf("prob_each_h: %d\t%f\n",i,prob_each_h[i]);}

	prob_each[0][0]=pressure_h2*3570*S; //H2 adsorption
	prob_each[0][1]=pressure_ch4*0.00874*S;//CH4 adsorption
	prob_each[0][2]=2.71*pow(10,13)*exp(-0.80/0.112);  //for H2 desorption
	prob_each[0][3]=2.71*pow(10,13)*exp(-0.58/0.112);  //for H attachment
	p[0][0]=-1;p[0][1]=1;// p for x even
	p[1][0]=-1;p[1][1]=0;
	p[2][0]=1;p[2][1]=0;
	q[0][0]=-1;q[0][1]=0;//q for x odd
	q[1][0]=1;q[1][1]=0;
	q[2][0]=1;q[2][1]=-1;
//	for(int i=row/2-G_row/2;i<row/2+G_row/2;i++)
//			for(int j=col/2-G_col/2;j<col/2+G_col/2;j++)
//				mesh[i][j][0]=9;
	for(int i=row/2-G_row/2-1;i<row/2+G_row/2;i++)
			for(int j=col/2-G_col/2-1;j<col/2+G_col/2;j++)
				mesh[i][j][0]=9;
	mesh[row/2-G_row/2-1][col/2-G_col/2-1][0]=100;
	mesh[row/2-G_row/2-1][col/2+G_col/2][0]=100;
	mesh[row/2+G_row/2][col/2-G_col/2-1][0]=100;
	mesh[row/2+G_row/2][col/2+G_col/2][0]=100;
}

int get_nb(int x,int y,int xy,int i){
	int x_temp,y_temp;
	if(xy==0){
		if(x%2==0)  x_temp=x+p[i][0];
		else        x_temp=x+q[i][0];
		if(x_temp>=row) return x_temp-=row;
		else if(x_temp<0) return x_temp+=row;
		else return x_temp;
	}
	else {
		if(x%2==0)  y_temp=y+p[i][1];
		else        y_temp=y+q[i][1];
		if(y_temp>=col) return y_temp-=col;
		else if(y_temp<0) return y_temp+=col;
		else return y_temp;
	}
}

void find_evt(){
	double prob_rand;
	double prob_sum=0;
	double prob_add=0;
	int drct,s_n;
	prob_sum+=prob_each[0][0];  //H2 adsorption
	prob_sum+=prob_each[0][1];  //CH4 adsorption
	prob_sum+=N_H*D_H*prob_each[0][2];  //H2 desorption
	//prob_sum+=N_H*D_H*prob_each[0][2];  //H2 desorption
//	prob_sum+=2*(G_row+G_col/2.0)*D_H*prob_each[0][3];  //H attach
	//prob_sum+=D_H*prob_each[0][3];  //H attach
	for(int i=0;i<9;i++)
		if(num_onsite[i]!=0) prob_sum+=num_onsite[i]*prob_onsite[i];
	for(int i=0;i<8;i++)
		if(!spec_list[i].empty())
			for(int j=0;j<spec_list[i].size();j++) prob_sum+=spec_list[i][j].prob_s;
	prob_rand=prob_sum*(double)rand()/RAND_MAX;
	prob_add+=prob_each[0][0];
	if(prob_add>prob_rand){
		num_evt[0]++;
		adsorption(1,1); //H2 adsorption
		//printf("adsorp h2\n");
		return;
	}
	prob_add+=prob_each[0][1];
	if(prob_add>=prob_rand){
		num_evt[2]++;
		adsorption(1,5); //CH4 adsorption
		//printf("adsorp ch4\n");
		return;
	}
	prob_add+=N_H*D_H*prob_each[0][2];
	if(prob_add>=prob_rand){
		num_evt[1]++;
		N_H-=2; //H2 desorption
		update_H();
		//printf("desorp h2\n");
		return;
	}
//	prob_add+=2*(G_row+G_col/2.0)*D_H*prob_each[0][3];
//	//prob_add+=D_H*prob_each[0][3];
//	if(prob_add>=prob_rand){
//		num_evt[21]++;
//		N_H--; //H attachment
//		update_H();
//		num_onsite[0]++;
//		attach_num[0]++;
//	//	printf("detach h\n");
//		return;
//	}

	for(int i=0;i<9;i++){ //detach
		if(num_onsite[i]!=0) {
			prob_add+=num_onsite[i]*prob_onsite[i];
			if(prob_add>=prob_rand) {
				num_onsite[i]--;
				num_det[i]++;
				if(i<4) {detach(i+1);num_evt[i+30]++;}
				else if(i<7)    {detach(i+2);num_evt[i+31]++;}
				else if(i==7)   {detach(3);num_detach_h[0]++;} //detachment of CH from H-saturated site
				else if(i==8)   {detach(6);num_detach_h[1]++;} //detachment of C2 from H-saturated site
		//		printf("detach cxhy\n");
				return;
			}
		}
	}
	for(int i=0;i<8;i++)
		if(!spec_list[i].empty()){
			for(int j=0;j<spec_list[i].size();j++){
				prob_add+=spec_list[i][j].prob_s;
				if(prob_add>=prob_rand){
					drct=spec_list[i][j].update();
					if(drct<3){
						if(spec_list[i][j].event[drct]==0){
							diffusion(i+1,j,drct);
							if(i>0) num_diff[i-1]++;
							//printf("diffuse cxhy");
						}
						else if(spec_list[i][j].event[drct]<=8){
							merge(i+1,j,drct);
							//printf("merge cxhy\n");
						}
						else if(spec_list[i][j].event[drct]==9){
							attach_num[i]++;
							if(i<5) num_onsite[i]++;
							else num_onsite[i-1]++;
							spec_del(i+1,j);
							//printf("attach cxhy\n");
							num_evt[i+21]++;
						}
						else if(spec_list[i][j].event[drct]==10){  //2017 added for attaching onto h-saturated sites
							if(i==2){
								attach_num_h[0]++;
								num_onsite[7]++;
								spec_del(i+1,j);
								//printf("attach cxhy\n");
								num_attach_h[0]++;
							}
							else if(i==5){
								attach_num_h[1]++;
								num_onsite[8]++;
								spec_del(i+1,j);
								//printf("attach cxhy\n");
								num_attach_h[1]++;
							}
							else err("error in new added");
						}
						else if(spec_list[i][j].event[drct]==100)
							return;
						else err("error in update1");
					}
					else  {
						spec_e(i+1,j,drct);
						if(evt_spec[i][drct-3]!=-1) num_evt[evt_spec[i][drct-3]-1]++;
					}
					return;
				}
			}
		}
	//err("prob");
}
void adsorption(int sp1,int sp2){
	int x,y,xn,yn;
	int flag=0;
	do{
		do{
			x=(int)(row*(double)rand()/RAND_MAX);
			y=(int)(col*(double)rand()/RAND_MAX);
		}while(mesh[x][y][0]!=0);
		for(int i=0;i<3;i++){    
			xn=get_nb(x,y,0,i);
			yn=get_nb(x,y,1,i);
			if(mesh[xn][yn][0]==0){
				flag=1;
				break;
			}
		}
	}while(flag==0);
//add 0913
	if(sp1==1){
		N_H++;
		update_H();
	}
	else creat(sp1,x,y);
	if(sp2>=0){
		if(sp2==1){
			N_H++;
			update_H();
		}
		else creat(sp2,xn,yn);
	}
}

void detach(int sp){
	int x,y,temp;
	temp=rand()%4;
	do{
		if(temp==0){
			x=row/2-G_row/2-2;
			y=col/2-G_col/2+(int)(G_col*(double)rand()/RAND_MAX);
		}
		else if(temp==1){
			x=row/2+G_row/2+1;
			y=col/2-G_col/2+(int)(G_col*(double)rand()/RAND_MAX);
		}
		else if(temp==2){
			x=row/2-G_row/2+2*(int)(G_row/2*(double)rand()/RAND_MAX);
			y=col/2-G_col/2-2;
		}
		else if(temp==3){
			x=row/2-G_row/2+1+2*(int)(G_row/2*(double)rand()/RAND_MAX);
			y=col/2+G_col/2+1;
		}
	}while(mesh[x][y][0]!=0);
	if(sp==1){
		N_H++;
		update_H();
	}
	else creat(sp,x,y);
}

void diffusion(int sp,int s_n,int drct){
	int x1,y1,x2,y2;
	x1=spec_list[sp-1][s_n].s_x;
	y1=spec_list[sp-1][s_n].s_y;
	x2=spec_list[sp-1][s_n].nb[drct][0];
	y2=spec_list[sp-1][s_n].nb[drct][1];
	if(mesh[x1][y1][0]!=sp||mesh[x2][y2][0]!=0){
		printf("in diffu,sp,x1,y1,x2,y2:  %d\t%d\t%d\t%d\t%d\n",sp,x1,y1,x2,y2);
		//if(mesh[x1][y1][0]!=sp) {printf("sp=%d,mesh=%d\n",sp,mesh[x1][y1][0]); err("diffusion1");}
		if(mesh[x1][y1][0]!=sp)  printf("dif1");
		if(mesh[x2][y2][0]!=0) printf("dif1");
			printf("drct=%d\tsp1=%d,mesh[0]=%d\n",drct,sp,mesh[x1][y1][0]);
			printf("event[drct]:%d",spec_list[sp-1][s_n].event[drct]);
			printf("sp2=%d,mesh2[1]=%d\n",mesh[x2][y2][0],mesh[x2][y2][1]);
			printf("x2=%d,y2[1]=%d\n",spec_list[sp-1][mesh[x2][y2][1]].s_x,spec_list[sp-1][mesh[x2][y2][1]].s_y);
			for(int i=0;i<3;i++){
				printf("evt[i]=%d\n",spec_list[sp-1][s_n].event[i]);
				update_evt_nb(x2,y2);
				printf("evt[i]=%d\n",spec_list[sp-1][s_n].event[i]);
				spec_list[sp-1][s_n].check_event();
				printf("evt[i]=%d\n",spec_list[sp-1][s_n].event[i]);
				printf("neigbour of 1: %d\t%d\n",get_nb(x1,y1,0,i),get_nb(x1,y1,1,i));
				printf("neigbour of 2: %d\t%d\n",get_nb(x2,y2,0,i),get_nb(x2,y2,1,i));
			}
			err("diffusion");
		//}
		err("diffusion");
	}
	mesh[x2][y2][0]=mesh[x1][y1][0];
	mesh[x2][y2][1]=mesh[x1][y1][1];
	mesh[x1][y1][0]=0;
	mesh[x1][y1][1]=0;
	spec_list[sp-1][s_n].s_x=x2;
	spec_list[sp-1][s_n].s_y=y2;
	spec_list[sp-1][s_n].check_event();
	update_evt_nb(x1,y1);
	update_evt_nb(x2,y2);
}

void merge(int sp1,int s_n,int drct){
	int sp2,sp_crt;
	int x1,y1,x2,y2;
	x1=spec_list[sp1-1][s_n].s_x;
	y1=spec_list[sp1-1][s_n].s_y;
	x2=spec_list[sp1-1][s_n].nb[drct][0];
	y2=spec_list[sp1-1][s_n].nb[drct][1];
	sp2=spec_list[sp1-1][s_n].event[drct];
	spec_del(sp1,s_n);
	spec_del(sp2,mesh[x2][y2][1]);
	sp_crt=sp_merge[sp1-1][sp2-1];
	if(evt_emerge[sp1-1][sp2-1]!=-1) num_evt[evt_emerge[sp1-1][sp2-1]-1]++;
	if(sp_crt>0) creat(sp_crt,x1,y1);
	else if(sp_crt==0) err("merge");
}

void spec_del(int sp,int s_n){
	int x,y;
	int spec_end;
	species_num[sp-1]--;
	x=spec_list[sp-1][s_n].s_x;
	y=spec_list[sp-1][s_n].s_y;
	if(sp!=mesh[x][y][0]) err("spec_del1");
	if(s_n!=mesh[x][y][1]) err("spec_del2");
	spec_end=spec_list[sp-1].size()-1;
	if(s_n!=spec_end){
	//	mesh[spec_list[sp-1][spec_end].s_x][spec_list[sp-1][spec_end].s_y][0]=sp;
		mesh[spec_list[sp-1][spec_end].s_x][spec_list[sp-1][spec_end].s_y][1]=s_n;
		swap(spec_list[sp-1][s_n],spec_list[sp-1][spec_end]);
	}
	spec_list[sp-1].pop_back();
	mesh[x][y][0]=0;
	mesh[x][y][1]=0;
	update_evt_nb(x,y);
}
void creat(int sp,int x,int y){
	species_num[sp-1]++;
	mesh[x][y][0]=sp;
	mesh[x][y][1]=spec_list[sp-1].size();
	spec_list[sp-1].push_back(species(sp,x,y));
	update_evt_nb(x,y);
}
void spec_e(int sp,int s_n,int drct){
	int sp1,sp2,xn,yn,flag=0;
	int x,y;
	x=spec_list[sp-1][s_n].s_x;
	y=spec_list[sp-1][s_n].s_y;
	sp1=sp_e[sp-1][drct-3][0];
	sp2=sp_e[sp-1][drct-3][1];
// revised 0925
	if(sp1==0)    err("spec_e");
	spec_del(sp,s_n);

	if(sp1==-2) {N_H++;update_H();}
	else if(sp1==-4) {N_H--;update_H();}
	else if(sp1>0) creat(sp1,x,y);
	
	if(sp2==-4) {N_H--;update_H();}
	else if(sp2>0){
		for(int i=0;i<3;i++){
			xn=get_nb(x,y,0,i);
			yn=get_nb(x,y,1,i);
			if(mesh[xn][yn][0]==0){
				creat(sp2,xn,yn);
				flag=1;
				break;
			}
		}
	}
	else if(sp2==-2) {N_H++;update_H();}
}

void update_evt_nb(int x,int y){
	int sp,s_n;
	int x_temp,y_temp;
	for(int i=0;i<3;i++){
		x_temp=get_nb(x,y,0,i);
		y_temp=get_nb(x,y,1,i);
		sp=mesh[x_temp][y_temp][0];
		if(sp>0&&sp<=8){
			s_n=mesh[x_temp][y_temp][1];
			spec_list[sp-1][s_n].check_event();
		}
	}
}

void statistic(int n,int N){
	int nn[8];
	for(int i=0;i<8;i++) nn[i]=spec_list[i].size();
	printf("%d+%d\t\tH\tC\tCH\tCH2\tCH3\tC2\tC2H\tC2H2\n",n,N);
	printf("N\t\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",species_num[0],species_num[1],species_num[2],species_num[3],species_num[4],species_num[5],species_num[6],species_num[7]);
	printf("N_nn\t\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",nn[0],nn[1],nn[2],nn[3],nn[4],nn[5],nn[6],nn[7]);
	printf("att\t\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",attach_num[0],attach_num[1],attach_num[2],attach_num[3],attach_num[4],attach_num[5],attach_num[6],attach_num[7]);
	printf("%d+%d\t\t\tH\tC\tCH\tCH2\tC2\tC2H\tC2H2\n",n,N);
	printf("det\t\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",num_det[0],num_det[1],num_det[2],num_det[3],num_det[4],num_det[5],num_det[6]);
	printf("onsite\t\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",num_onsite[0],num_onsite[1],num_onsite[2],num_onsite[3],num_onsite[4],num_onsite[5],num_onsite[6]);
	printf("num_H:\t%f\n",N_H);
}

int get_sp(int x,int y){
	return mesh[x][y][0];
}

void update_H(){
	D_H=N_H/S;
	for(int i=1;i<7;i++)
		if(!spec_list[i].empty())
			for(int j=0;j<spec_list[i].size();j++)
				spec_list[i][j].upd_H();
}

