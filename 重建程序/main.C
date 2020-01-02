#include <iostream>

#include <vector>

#include <fstream>

#include <sstream>

#include <math.h>

#include <stdlib.h>

#include <TMath.h>
#include <TGraph2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF2.h>
#include <TH1.h>
#include <Math/Functor.h>
#include <TPolyLine3D.h>
 #include <Math/Vector3D.h>
#include <Fit/Fitter.h>
//#include <Eigen/Dense>
#include <cassert>

using namespace ROOT::Math;

using namespace std;



int number[40]; //记录筛选出来的事件编号  265053

int x[40][48],y01[40][40],y2[40][40],a1[40][40],a2[40][40],a3[40][48];  //读取数据

double x_1[120],y_1[120],z_1[120],a_1[120]; //记录重建点的坐标及adc值

//double x_2[50],y_2[50],z_2[50];

double velocity[40]={10.143675,9.94548,10.1305825,10.0936525,9.85208,9.881315,9.9555725,9.99623,10.263805,10.3858125,10.2523,10.2861,10.073725,10.3249,10.345685,10.21335,10.42214,10.363475,10.264895,10.40365,10.299835,10.258625,10.142685,10.127855,9.927885,10.121475,10.175885,10.1780525,10.052285,9.909235,10.14835,9.9660725,10.017645,9.9072075,9.9061125,9.769795,9.919345,9.8669525,9.723265,9.8208025};

double X[40][120],Y[40][120],Z[40][120],A[40][120];

int N[40];

void line(double t, const double *p, double &x, double &y, double &z)   //定义直线方程，四个待定系数
//void line(double t, double *p, double &x, double &y, double &z)
{      
   x = p[0] + p[1]*t;
   y = p[2] + p[3]*t;
   z = t;
}

bool first = true;

struct SumDistance2 {                                     //算所有散点到某直线距离的和（我自己也不太懂，直接复制过来的）
   // the TGraph is a data member of the object
   TGraph2D *fGraph;

   SumDistance2(TGraph2D *g) : fGraph(g) {}

   // calculate distance line-point
   double distance2(double x,double y,double z, const double *p) {
   //double distance2(double x,double y,double z, double *p) {
      // distance line point is D= | (xp-x0) cross  ux |
      // where ux is direction of line and x0 is a point in the line (like t = 0)
      XYZVector xp(x,y,z);
      XYZVector x0(p[0], p[2], 0. );
      XYZVector x1(p[0] + p[1], p[2] + p[3], 1. );
      XYZVector u = (x1-x0).Unit();
      double d2 = ((xp-x0).Cross(u)).Mag2();
      return d2;
   }

   // implementation of the function to be minimized
   double operator() (const double *par) {
   //double operator() (double *par) {
      assert(fGraph != 0);
      double * x = fGraph->GetX();
      double * y = fGraph->GetY();
      double * z = fGraph->GetZ();
      int npoints = fGraph->GetN();
      double sum = 0;
      for (int i  = 0; i < npoints; ++i) {
         double d = distance2(x[i],y[i],z[i],par);
         sum += d;
      }
      if (first) {
         std::cout << "Total Initial distance square = " << sum << std::endl;
      }
      first = false;
      return sum;
   }

};


Int_t line3Dfit(double *x,double *y,double *z,int p,double *par,double *ps)  //散点3维拟合
{
   gStyle->SetOptStat(0);
   gStyle->SetOptFit();


   //double e = 0.1;
   Int_t nd = 10000;


   // double xmin = 0; double ymin = 0;
   // double xmax = 10; double ymax = 10;

  TGraph2D * gr = new TGraph2D();

for(int N=0;N<p;N++)
{
	gr->SetPoint(N,x[N],y[N],z[N]);
	//cout<<N<<"\t"<<x[N]<<"\t"<<y[N]<<"\t"<<z[N]<<endl;
}
   
   ROOT::Fit::Fitter  fitter;


   // make the functor objet
   SumDistance2 sdist(gr);            //算刚刚设置点到已知直线的距离和 
   ROOT::Math::Functor fcn(sdist,4);
   // set the function and the initial parameter values
   double pStart[4] = {ps[0],ps[1],ps[2],ps[3]};
   fitter.SetFCN(fcn,pStart);
   // set step sizes different than default ones (0.3 times parameter values)
   for (int i = 0; i < 4; ++i) fitter.Config().ParSettings(i).SetStepSize(0.01);

   bool ok = fitter.FitFCN();
   if (!ok) {
      Error("line3Dfit","Line3D Fit failed");
      return 1;
   }

   const ROOT::Fit::FitResult & result = fitter.Result();

   std::cout << "Total final distance square " << result.MinFcnValue() << std::endl;
   result.Print(std::cout);


   //gr->Draw("p0");

   const double * parFit = result.GetParams();
   //double * parFit = result.GetParams();

    parFit=result.GetParams();
for(int a=0;a<4;a++)
{
	*(par+a)=*(parFit+a);
}
    
	return 0;
}  

void DrawFitline (double *parFit)                       //画出拟合后的线
{
 /*  int n = 1000;
   double t0,dt;            //确定画图时t的范围，t0增大到dt
   if(-parFit[0]/parFit[1]<(100-parFit[0])/parFit[1])    //通过x的范围确定t的范围
   {
	t0=-parFit[0]/parFit[1];
	dt=(100-parFit[0])/parFit[1];
   }
   if(-parFit[0]/parFit[1]>=(100-parFit[0])/parFit[1]) 
   {
	t0=(100-parFit[0])/parFit[1];
	dt=-parFit[0]/parFit[1];
   }
   TPolyLine3D *l = new TPolyLine3D(1000);
   for (int i = 0; i <n;++i) {
      double t = t0+ dt*i/n;
      double x,y,z;
      line(t,parFit,x,y,z);
      l->SetPoint(i,x,y,z);  
   }   */
   double t0,dt;
   t0=-parFit[0]/parFit[1];
   dt=(100-parFit[0])/parFit[1];
   TPolyLine3D *l = new TPolyLine3D();
   double x,y,z;
   line(t0,parFit,x,z,y);//line(t0,parFit,x,y,z);
   l->SetPoint(0,x,y,z); //l->SetPoint(0,x,y,z);
   line(dt,parFit,x,z,y);//line(dt,parFit,x,y,z);
   l->SetPoint(1,x,y,z); //l->SetPoint(1,x,y,z);


   l->SetLineColor(kRed);
   l->Draw("same");
}

void DrawFitline1 (double *parFit)                       //画出拟合后的线
{   double t0,dt;
   t0=0;
   dt=25.92;
   TPolyLine3D *l = new TPolyLine3D();
   double x,y,z;
   line(t0,parFit,x,z,y);//line(t0,parFit,x,y,z);
   l->SetPoint(0,x,y,z); //l->SetPoint(0,x,y,z);
   line(dt,parFit,x,z,y);//line(dt,parFit,x,y,z);
   l->SetPoint(1,x,y,z); //l->SetPoint(1,x,y,z);


   l->SetLineColor(kRed);
   l->SetLineWidth(3);
   l->Draw("same");
}
void DrawFitline2 (double *parFit)                       //画出拟合后的线
{   double t0,dt;
   t0=0;
   dt=25.92;
   TPolyLine3D *l = new TPolyLine3D();
   double x,y,z;
   line(t0,parFit,x,z,y);//line(t0,parFit,x,y,z);
   l->SetPoint(0,x,y,z); //l->SetPoint(0,x,y,z);
   line(dt,parFit,x,z,y);//line(dt,parFit,x,y,z);
   l->SetPoint(1,x,y,z); //l->SetPoint(1,x,y,z);


   l->SetLineColor(kBlue);
   l->SetLineWidth(3);
   l->Draw("same");
}



//寻找碰撞点坐标（学长写的，没有在用

int find_vertex(double *x_1,int *y_1, double *z_1, double *a_1, int cnt){

  for(int i=cnt-1;i>0;--i){

    //cout<<a_1[i-1]-a_1[i]<<" ("<<y_1[i]<<") ";

    if((a_1[i-1]-a_1[i]) > 7500) { return i; }

  }

  return -1;

}





void reconstruct(int number, int *x, int *y1, int *a1,int *a3,int K){

  //cout<<number<<endl;



  int cnt=0;  // 第几个重建的点

  int xmax=-9999,xmin=9999,ymax=-9999,ymin=9999;

  int index_xmax=0, index_xmin=0, index_ymax=0, index_ymin=0;

  int error=0;  //允许的同一个点的时间误差

  

  //找到x中的最大、最小值，并记录索引

  for(int k=0; k<48; ++k){

    if(x[k]>xmax){

      xmax=x[k];

      index_xmax=k;

    }

    if(x[k]<xmin && x[k]>0){

      xmin=x[k];

      index_xmin=k;

    }

  }

  

  //找到y中的最大、最小值，并记录索引

  for(int j=0; j<40; ++j){

    if(y1[j]>ymax){

      ymax=y1[j];

      index_ymax=j;

    }

    if(y1[j]<ymin && y1[j]>0){

      ymin=y1[j];

      index_ymin=j;

    }

  }

  

   

  //cout<<"y->x"<<endl;

  //在x中找到与各个y对应的时间最接近的两个时间，记录这两个x的索引

  for(int j=0; j<40; ++j){

    if(y1[j]>=xmin-error && y1[j]<=xmax+error){

      int up=9999, down=9999;

      int index_up=index_xmax, index_down=index_xmin; //将两个索引初始值设为最大最小值的索引

      

      for(int k=0; k<48; ++k){

        if(x[k]>0){

          if(y1[j]==x[k]){

            index_up=index_down=k;

            up=down=0;

          }

          if((x[k]>y1[j]) && ((x[k]-y1[j])<up)){

            up=x[k]-y1[j];

            index_up=k;

          }

          else if((x[k]<y1[j]) && ((y1[j]-x[k])<down)){

            down=y1[j]-x[k];

            index_down=k;

          }

        }

      }

      

      y_1[cnt]=j;

      if(index_down == index_up){

        x_1[cnt]=index_down;

      }

      else{

        x_1[cnt]=index_down+(index_up-index_down)*(y1[j]-x[index_down])/(double)(x[index_up]-x[index_down]);

      }

      z_1[cnt]=y1[j];

      a_1[cnt]=a1[j];

      x_1[cnt]=x_1[cnt]*2.08333;

      y_1[cnt]=y_1[cnt]*2.5;

      //z_1[cnt]=z_1[cnt]/100*velocity[j];
	z_1[cnt]=z_1[cnt]/10;
       X[K][N[K]]=x_1[cnt];  //就是把数据重新整合到一个X，Y，Z，A的数组里面

	  Y[K][N[K]]=y_1[cnt];

	  Z[K][N[K]]=z_1[cnt];
	  A[K][N[K]]=a_1[cnt];

	  N[K]++;

      //cout<<x_1[cnt]<<"\t"<<y_1[cnt]<<"\t"<<z_1[cnt]<<endl;

      ++cnt;

    }

  }

 

  //cout<<"x->y"<<endl;

  //在y中找到与各个x对应的时间最接近的两个时间，记录这两个y的索引

 for(int k=0; k<48; ++k){

    if(x[k]>=ymin-error && x[k]<=ymax+error){

      int up=9999,down=9999;

      int index_up=index_ymax, index_down=index_ymin;

      

      for(int j=0; j<40; ++j){

        if(y1[j]>0){

          if(y1[j]==x[k]){

            index_up=index_down=j;

            up=down=0;

          }

          else if((y1[j]>x[k]) && ((y1[j]-x[k])<up)){

            up=y1[j]-x[k];

            index_up=j;

          }

          else if((y1[j]<x[k]) && ((x[k]-y1[j])<down)){

            down=x[k]-y1[j];

            index_down=j;

          }

        }

      }

      

      x_1[cnt]=k;

      if(index_down == index_up){

        y_1[cnt]=index_down;

      }

      else{

        y_1[cnt]=index_down+(index_up-index_down)*(x[k]-y1[index_down])/(double)(y1[index_up]-y1[index_down]);

      }

      z_1[cnt]=x[k];

      a_1[cnt]=a3[k];

      x_1[cnt]=x_1[cnt]*2.08333;

      //z_1[cnt]=z_1[cnt]/100*velocity[int(floor(y_1[cnt]))];//y_1[cnt]可能是小数
	z_1[cnt]=z_1[cnt]/10;
      y_1[cnt]=y_1[cnt]*2.5;

      X[K][N[K]]=x_1[cnt];

	  Y[K][N[K]]=y_1[cnt];

	  Z[K][N[K]]=z_1[cnt];
	  A[K][N[K]]=a_1[cnt];
	  N[K]++;

      //cout<<x_1[cnt]<<"\t"<<y_1[cnt]<<"\t"<<z_1[cnt]<<endl;

      ++cnt;

    }

  } 

 

}




double D(double x1,double y1,double z1,double x2,double y2,double z2) //计算两点间的距离
{
	double d=0;
	d=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
	return d;
}


void Sort (double *X1, double *Y1, double *Z1,double *A1) //  把所有数据按y值从大到小排列（粒子是在y减小的方向入射）
{
        double k = 0;
	for (int i = 0; i< 120; i++)
	{
		for (int j = i;j<120; j++)
		{
			if(Y1[i]<Y1[j])
			{
				
				k = Y1[j];
				Y1[j] = Y1[i];
				Y1[i] = k;
				k = X1[j];
				X1[j] = X1[i];
				X1[i] = k;
				k = Z1[j];
				Z1[j] = Z1[i];
				Z1[i] = k;
				k = A1[j];
				A1[j] = A1[i];
				A1[i] = k;
			}
		}
	}
	
}

void separate(int m,double *M,double *X1, double *Y1, double *Z1,double *A1,double *X2,double *Y2,double *Z2,double *X3,double *Y3,double *Z3,int *pp,int *qq)//初步找寻径迹分岔处
{
	int i;
	int p=0,q=0;
	for(i=0;i<m;i++)
	{
		if(Y1[i]<=M[1])
		break;
	}
	for(int j=i;j<m;j++)
		{
			if(X1[j]<=M[0]) //左分岔
			  {
				X2[p]=X1[j];
				Y2[p]=Y1[j];
				Z2[p]=Z1[j];
				p++;
			  }
			else{                          //右分岔
				X3[q]=X1[j];
				Y3[q]=Y1[j];
				Z3[q]=Z1[j];
				q++;
			    }
		}
	*pp = p;
	*qq = q;
}


int findpoint(int m,double *X1, double *Y1, double *Z1,double *A1,double *X2,double *Y2,double *Z2,double *X3,double *Y3,double *Z3,int *pp,int *qq)//初步找寻径迹分岔处
{       int i=0;
        int p = 0;
	int q = 0;
	int en=10000;
	int e=15;//定义的相邻点之间最小分开距离，单位mm
	for(i=0;i<m;i++)
	{
			//if(fabs(*(X1+i)-*(X1+i+1))>=e||fabs(*(X1+i)-*(X1+i+2))>=e||fabs(*(X1+i)-*(X1+i+3))>=e)
			if((A1[i]>=en||A1[i+1]>=en||A1[i+2]>=en)&&(D(X1[i],Y1[i],Z1[i],X1[i+1],Y1[i+1],Z1[i+1])>=e||D(X1[i],Y1[i],Z1[i],X1[i+2],Y1[i+2],Z1[i+2])>=e))
			
			 
			break;
	}
	for(int j=i;j<m;j++)
		{
			if(X1[j]>=(X1[i]+X1[i+1]+X1[i+2])/3) //左分岔
			  {
				X2[p]=X1[j];
				Y2[p]=Y1[j];
				Z2[p]=Z1[j];
				p++;
			  }
			else{                          //右分岔
				X3[q]=X1[j];
				Y3[q]=Y1[j];
				Z3[q]=Z1[j];
				q++;
			    }
		}
	*pp = p;
	*qq = q;
	return i;
			
}

double footpoint1(double *par1,double *par2,double *M)  //求两条直线公垂线中点坐标
{
	double t1,t2,x1,y1,z1,x2,y2,z2;
	double a,b,c,a1,b1,c1;
	double d;
	a=par1[1]*par1[1]+par1[3]*par1[3]+1;
	b=-(par1[1]*par2[1]+par1[3]*par2[3]+1);
	c=par1[1]*(par1[0]-par2[0])+par1[3]*(par1[2]-par2[2]);
	a1=par1[1]*par2[1]+par1[3]*par2[3]+1;
	b1=-(par2[1]*par2[1]+par2[3]*par2[3]+1);
	c1=par2[1]*(par1[0]-par2[0])+par2[3]*(par1[2]-par2[2]);
	t1=(c1*b-c*b1)/(a*b1-a1*b);
	t2=(c*a1-c1*a)/(a*b1-a1*b);
	line(t1,par1,x1,y1,z1);
	line(t2,par2,x2,y2,z2);
	M[0]=(x1+x2)/2;
	M[2]=(y1+y2)/2;
	M[1]=(z1+z2)/2;
	d=D(x1,y1,z1,x2,y2,z2);
	return d;
}
/*void footpoint2(double *par,double *x,double *y,double *z,int n,)//求点到直线距离最近的点,未完待续。。。
{
	double t,d,m;
	double x0,y0,z0;
	d=100;
	for(int i=0;i<n;i++)
	{
		t=(z[i]-par[1]*(par[0]-x[i])-par[3]*(par[2]-y[i]))/(par[1]*par[1]+par[3]*par[3]+1);
		line(t,par,x0,y0,z0);
		if(d>=D(x[i],y[i],z[i],x0,y0,z0))
		{
			d=D(x[i],y[i],z[i],x0,y0,z0);
			m=i;
		
	
*/


double angle(double *p1,double *p2)
{
	double angle=acos((p1[1]*p2[1]+p1[3]*p2[3]+1)/(D(p1[1],p1[3],1,0,0,0)*D(p2[1],p2[3],1,0,0,0)));
	angle=angle*180/3.1415926;
	return angle;
}

int main(int nstart){



for(int s=0;s<40;s++)

 {N[s]=0;}

  int index_branch=-1;  //第几个筛选出来的事件

  int switch_leaf=0;  //事件中的第几类数据

  ifstream infile;

  infile.open("/home/jeff/TPC/data_new.txt"); //改成你的路径

  //string strTmp;

  long iTmp; //读到的数

  for(int i=0; i<6682;++i){

    infile>>iTmp;

  

/*  while(getline(infile,strTmp,' ')){

    istringstream iss(strTmp);

    iss >> iTmp;*/

    if(iTmp<0) { iTmp=-1; }

    int jaziel = i%257; //每个事件包含209个数，jaziel记录该数在事件中的序号

    if(jaziel == 0){  //第一个数为事件编号

      ++index_branch;

      switch_leaf=0;

      number[index_branch]=iTmp;

      //cout<<i<<"  "<<iTmp<<endl;

    }

    

    else if(jaziel<=40)switch_leaf=1; //第2~41个数为y1_tdc

    else if(jaziel<=80)switch_leaf=2; //第42~81个数为y2_tdc

    else if(jaziel<=120)switch_leaf=3;//第82~121个数为y1_adc

    else if(jaziel<=160)switch_leaf=4;//第122~161个数为y2_adc

    else if(jaziel<=208)switch_leaf=5;//第162~209个数为x_tdc

    else switch_leaf=6;               //第210~257个数为x_adc

    switch(switch_leaf){

      case 0:

        break;

      case 1:

        y01[index_branch][jaziel-1]=iTmp;

        break;

      case 2:

        y2[index_branch][jaziel-41]=iTmp;
	break;

      case 3:

        a1[index_branch][jaziel-81]=iTmp;

        break;

      case 4:

        a2[index_branch][jaziel-121]=iTmp;

        break;

      case 5:

        x[index_branch][jaziel-161]=iTmp;

        break;

      case 6:

        a3[index_branch][jaziel-209]=iTmp;

        break;

      default:

        cout<<jaziel<<"  error!!"<<endl;

    }


  }

  //TGraph2D *g;

 /* for(int a=0;a<20;a++)

 {

  g[a]=new TGraph2D();

 }*/

  TCanvas *c1=new TCanvas("c1","",600,600);
  TCanvas *c2=new TCanvas("c2","",600,600);
  TGraph2D *g=new TGraph2D();
  TCanvas *c3=new TCanvas("c3","",600,600);

  c1->SetMargin(0.2,0.2,0.1,0.1);
  c3->Divide(2,2);
  	TH2F *hadcx = new TH2F("hadcx","adc vs X ch",48,-0.5,47.5,200,50,46050);
	hadcx->SetStats(false);
	TH2F *hadcy1 = new TH2F("hadcy1","adc vs Y1 ch",40,-0.5,39.5,200,50,46050);
	hadcy1->SetStats(false);
	//TH2F *hadcy2 = new TH2F("hadcy2","adc vs Y2 ch",40,39.5,79.5,200,1,101);
	//hadcy2->SetStats(false);
	//TH2F *htdcx = new TH2F("htdcx","tdc vs X ch",48,-0.5,47.5,200,-000,2000);
	TH2F *htdcx = new TH2F("htdcx","tdc vs X ch",48,0,100,200,-000,140);
	htdcx->SetStats(false);
   	htdcx->SetMarkerStyle(5);
   	htdcx->SetMarkerSize(0.7);
   	htdcx->SetMarkerColor(9);
	//TH2F *htdcy1 = new TH2F("htdcy1","tdc vs Y1 ch",40,-0.5,39.5,200,-000,100);
	TH2F *htdcy1 = new TH2F("htdcy1","tdc vs Y1 ch",40,0,100,200,-000,140);
	htdcy1->SetStats(false);
   	htdcy1->SetMarkerStyle(5);
   	htdcy1->SetMarkerSize(0.7);
   	htdcy1->SetMarkerColor(9);
	/*TH2F *htdcy2 = new TH2F("htdcy2","tdc vs Y2 ch",40,39.5,79.5,200,-000,2000);
	htdcy2->SetStats(false);
   	htdcy2->SetMarkerStyle(5);
   	htdcy2->SetMarkerSize(0.7);
   	htdcy2->SetMarkerColor(9);*/

        

  for(int k=nstart; k<=index_branch; ++k)

  {

   // getchar();
    double X2[120],Y2[120],Z2[120],X3[120],Y3[120],Z3[120];

    int p=0,q=0,n=0;

    reconstruct(number[k],x[k],y01[k],a1[k],a3[k],k);

    reconstruct(number[k],x[k],y2[k],a2[k],a3[k],k);
 
    Sort(X[k],Y[k],Z[k],A[k]);     //把点按照y从大到小的顺序排列，即粒子前进的方向

    
    n=findpoint(N[k],X[k],Y[k],Z[k],A[k],X2,Y2,Z2,X3,Y3,Z3,&p,&q); //初步寻找径迹分岔处，返回找的分岔点序号

    cout<<"constructing "<<k<<"........"<<endl;

    cout<<number[k]<<endl;

    g->SetTitle(";x/mm;y/mm;z/mm");

    //g->SetMaximum(1000);



    g->SetPoint(0,0,0,0);

    g->SetPoint(1,100,100,140);

    //g->Clear();   

//line3Dfit(double *x,double *y,double *z,int n)

/*for(int j=0;j<q;j++)     //分开后左支（右支为q）的点
{
	//cout<<X2[j]<<"\t"<<Y2[j]<<"\t"<<Z2[j]<<"\t"<<endl;
	//g->SetPoint(j+2,X2[j],Y2[j],Z2[j]);
	g->SetPoint(j+2,X3[j],Y3[j],Z3[j]);
} */

/*for(int j=n;j<N[k];j++)   //分开后的两支的点
{
	g->SetPoint(j-n+2,X[k][j],Y[k][j],Z[k][j]);
		
}  */
 for(int m=0;m<N[k];m++)  //第k组所有的点

    {
         
       cout<<X[k][m]<<"\t"<<Y[k][m]<<"\t"<<Z[k][m]<<"\t"<<A[k][m]<<endl;

       g->SetPoint(m+2,X[k][m],Y[k][m],Z[k][m]);

		

    } 

  	c1->cd();

	g->Draw(); 
                      
	double par1[4],par2[4],M[3],ps1[4],ps2[4];
	
	//ps数组为par拟合初始值，这里取数组中始末两点的连线计算初始值
	ps1[1]=(X2[p-1]-X2[0])/(Y2[p-1]-Y2[0]);
	ps1[0]=X2[p-1]-ps1[1]*Y2[p-1];
	ps1[3]=(Z2[p-1]-Z2[0])/(Y2[p-1]-Y2[0]);
	ps1[2]=Z2[p-1]-ps1[3]*Y2[p-1];

	ps2[1]=(X3[q-1]-X3[0])/(Y3[q-1]-Y3[0]);
	ps2[0]=X3[q-1]-ps2[1]*Y3[q-1];
	ps2[3]=(Z3[q-1]-Z3[0])/(Y3[q-1]-Y3[0]);
	ps2[2]=Z3[q-1]-ps2[3]*Y3[q-1];
                      
	line3Dfit(X2,Z2,Y2,p,par1,ps1); //拟合左半部分散点
	
	//cout<<"lalala"<<endl;
	//DrawFitline (par1);          //画出左半部分的拟合线
	
	
	
	line3Dfit(X3,Z3,Y3,q,par2,ps2); //拟合右半部分散点                     
	//cout<<"bababa"<<endl;
	//DrawFitline (par2);          //画出右半部分的拟合线
	
	double d;
	d=footpoint1(par1,par2,M);      //求出拟合线公垂线的中点，存在M[3]中


	separate(N[k],M,X[k],Y[k],Z[k],A[k],X2,Y2,Z2,X3,Y3,Z3,&p,&q);
        
	ps1[1]=(X2[p-1]-X2[0])/(Y2[p-1]-Y2[0]);
	ps1[0]=X2[p-1]-ps1[1]*Y2[p-1];
	ps1[3]=(Z2[p-1]-Z2[0])/(Y2[p-1]-Y2[0]);
	ps1[2]=Z2[p-1]-ps1[3]*Y2[p-1];

	ps2[1]=(X3[q-1]-X3[0])/(Y3[q-1]-Y3[0]);
	ps2[0]=X3[q-1]-ps2[1]*Y3[q-1];
	ps2[3]=(Z3[q-1]-Z3[0])/(Y3[q-1]-Y3[0]);
	ps2[2]=Z3[q-1]-ps2[3]*Y3[q-1];
	
	line3Dfit(X2,Z2,Y2,p,par1,ps1);//第二次拟合
	DrawFitline1 (par1);
	line3Dfit(X3,Z3,Y3,q,par2,ps2);
	DrawFitline2 (par2);

	double a=angle(par1,par2);
	cout<<"angle="<<a<<endl;
	cout<<"交点坐标为("<<M[0]<<","<<M[1]<<","<<M[2]<<")"<<endl;
	cout<<"两直线之间距离为"<<d<<endl;

	TGraph2D *g1=new TGraph2D();
	g1->SetPoint(0,0,0,0);

    	g1->SetPoint(1,100,100,140);
	for(int j=0;j<p;j++)     //分开后左支（右支为q）的点

{
	g1->SetPoint(j+2,X2[j],Y2[j],Z2[j]);
}
	for(int j=0;j<q;j++)     //分开后左支（右支为q）的点

{
	g1->SetPoint(j+p+2,X3[j],Y3[j],Z3[j]);
}

	//g->SetPoint(100,M[0],M[1],M[2]);
	
	c2->cd();	
	g1->Draw();

 	for(int l=0;l<40;l++){

					hadcy1->Fill(l,a1[k][l]);
					//htdcy1->Fill(l,tdc_y[l]);
					hadcy1->Fill(l,a2[k][l]);
					htdcy1->Fill(float(l+1)*2.5,float(y01[k][l])/10);
					htdcy1->Fill(float(l+1)*2.5,float(y2[k][l])/10);
				}
	for(int j=0;j<48;j++)
	{
	hadcx->Fill(j,a3[k][j]);
	htdcx->Fill(float(j+1)*2.08333,float(x[k][j])/10);
	}
	gStyle->SetPalette(100);
        c3->cd(1);
	hadcx->Draw("colz");
        c3->cd(2);
	hadcy1->Draw("colz");  
        c3->cd(3);
	htdcx->Draw("colz");  
        c3->cd(4);
	htdcy1->Draw("colz");                                                          
	
    //c1->cd();

                                      //g->Draw();                                                             
//第八组 第十组
   // c1->Update();

                //0 90.2221;1 85.8157;3 88
 
       				//13,16,17,19（缺少中间，对的很好）,21,22,23 ,24            	 
 break; 
                                                                                                                              

  }

/*ofstream outfile;                  //把事件号和所有点的坐标输出到data.txt文件中

 outfile.open("data.txt");

for(int b=0;b<=index_branch; ++b)

{



  outfile<<"\n"<<number[b]<<"\n"<<endl;

  for(int e=0;e<N[b];e++)

  {

    outfile<<X[b][e]<<"\t"<<Y[b][e]<<"\t"<<Z[b][e]<<"\t"<<endl;

  }

}

outfile.close();                   */                                    



}

// 2 5 12
