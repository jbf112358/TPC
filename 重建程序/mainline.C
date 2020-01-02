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
#include <cassert>
#include "SetHstyle.h"
#define event 13535    //事件数目 
using namespace ROOT::Math;

using namespace std;


 
int number[event]; //记录筛选出来的事件编号
int Fail;
//int x[event][48],y01[event][40],y2[event][40],a1[event][40],a2[event][40],a3[event][48];  //读取数据
int x[event][48],y01[event][40],y2[event][40],a1[event][40],a2[event][40],a3[event][48];
double x_1[120],y_1[120],z_1[120],a_1[120]; //记录重建点的坐标及adc值

//double x_2[50],y_2[50],z_2[50];

double velocity[80]={9.91872,9.98426,10.004,9.9719,10.0258,10.0368,10.0297,9.97305,10.0361,9.97746,10.0105,9.95898,10.0505,10.0446,10.0392,9.95529,10.0149,9.95702,10.0023,10.0465,10.0438,9.96387,9.99135,10.0069,10.0014,9.9962,9.95366,9.98906,10.0075,9.98603,10.011,10.0338,10.0009,9.96944,10.0291,9.98447,10.0335,10.0123,10.0286,9.91186,9.9347,9.96115,9.95841,10.0022,10.0691,10.0193,10.036,10.0237,10.0008,10.0095,10.0049,9.98606,10.0324,9.98899,10.001,9.97822,9.96624,9.96806,9.98085,9.9923,9.94606,10.0356,10.044,10.0292,10.0787,10.0249,10.004,10.0029,10.0101,10.0257,9.97714,9.97251,9.99361,10.0156,9.98248,10.0367,9.96419,9.97412,9.9632,10.0006};
double resmean[80]={5.73645,1.1032,-0.279141,1.97273,-1.80419,-2.56439,-2.07343,1.89151,-2.52102,1.58151,-0.735298,2.88295,-3.51541,-3.10488,-2.73035,3.14367,-1.04274,3.02169,-0.164296,-3.24289,-3.04966,2.53847,0.605954,-0.479503,-0.0980745,0.266408,3.25892,0.766844,-0.527783,0.978944,-0.768808,-2.35656,-0.066097,2.14605,-2.02774,1.0891,-2.33735,-0.859667,-1.99446,6.22444,4.60079,2.73008,2.92344,-0.156324,-4.80287,-1.34627,-2.5119,-1.65835,-0.056141,-0.663348,-0.34386,0.977076,-2.26105,0.771831,-0.0695591,1.52771,2.37129,2.2432,1.3434,0.539295,3.79643,-2.4854,-3.06826,-2.037,-5.46908,-1.74115,-0.282715,-0.205019,-0.707683,-1.79614,1.60397,1.92935,0.447425,-1.08846,1.22858,-2.55658,2.51576,1.81646,2.58522,-0.044055};
double X[event][120],Y[event][120],Z[event][120],A[event][120],a[event],d[event],leth[event];
double mc[event];//每个事件的y-x误差
int N[event];
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

  /* std::cout << "Total final distance square " << result.MinFcnValue() << std::endl;
   result.Print(std::cout);*/


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

void DrawFitline1 (double *parFit)                       //画出拟合后前半段的线
{
   double t0,dt;
   t0=100;
   dt=50;//50
   TPolyLine3D *l = new TPolyLine3D();
   double x,y,z;
   line(t0,parFit,x,z,y);//line(t0,parFit,x,y,z);
   l->SetPoint(0,x,y,z); //l->SetPoint(0,x,y,z);
   line(dt,parFit,x,z,y);//line(dt,parFit,x,y,z);
   l->SetPoint(1,x,y,z); //l->SetPoint(1,x,y,z);

   l->SetLineWidth(3);
   l->SetLineColor(kRed);
   l->Draw("same");
}
void DrawFitline2 (double *parFit)                       //画出拟合后后半段的线
{
   double t0,dt;
   t0=50;//50
   dt=0;
   TPolyLine3D *l = new TPolyLine3D();
   double x,y,z;
   line(t0,parFit,x,z,y);//line(t0,parFit,x,y,z);
   l->SetPoint(0,x,y,z); //l->SetPoint(0,x,y,z);
   line(dt,parFit,x,z,y);//line(dt,parFit,x,y,z);
   l->SetPoint(1,x,y,z); //l->SetPoint(1,x,y,z);

   l->SetLineWidth(3);
   l->SetLineColor(kBlue);
   l->Draw("same");
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
double footpoint(double *par1,double *par2,double *M)  //求两条直线公垂线中点坐标
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

void findpoint(int m,double *X1, double *Y1, double *Z1,double *A1,double *X2,double *Y2,double *Z2,double *X3,double *Y3,double *Z3,int *pp,int *qq)//按平面把轨迹截断，同时去掉adc<1000的，分开的两端存在X2,X3中
{   
        int p = 0;
	int q = 0;
	int e=50;//定义的分开直线的面，单位mm
	
	for(int j=0;j<m;j++)
		{
			if(Y1[j]>e) //前半段
			  {
				X2[p]=X1[j];
				Y2[p]=Y1[j];
				Z2[p]=Z1[j];
				p++;
			  }
			else{                          
				X3[q]=X1[j];
				Y3[q]=Y1[j];
				Z3[q]=Z1[j];
				q++;
			    }
		}
	*pp = p;
	*qq = q;
			
}
void findpoint1(int m,double *X1, double *Y1, double *Z1,double *A1,double *X2,double *Y2,double *Z2,double *X3,double *Y3,double *Z3,int *pp,int *qq)//把轨迹按隔一个点分开,分开的两部分存在X2，X3中
{
        int p = 0;
	int q = 0;
	for(int j=0;j<m;j++)
	{
		if(j%2==1)
		{
			X2[p]=X1[j];
			Y2[p]=Y1[j];
			Z2[p]=Z1[j];
			p++;
		}
		else
		{
			X3[q]=X1[j];
			Y3[q]=Y1[j];
			Z3[q]=Z1[j];
			q++;
		}
	}
	*pp = p;
	*qq = q;
}
double angle(double *p1,double *p2)
{
	double angle=acos((p1[1]*p2[1]+p1[3]*p2[3]+1)/(D(p1[1],p1[3],1,0,0,0)*D(p2[1],p2[3],1,0,0,0)));
	angle=angle*180/3.1415926;
	return angle;
}

double length(double *p1,double *p2)
{
	double leth,t0=50;
   	TPolyLine3D *l = new TPolyLine3D();
   	double x1,y1,z1,x2,y2,z2;
   	line(t0,p1,x1,z1,y1);//line(t0,parFit,x,y,z);
	line(t0,p2,x2,z2,y2);
	leth=D(x1,y1,z1,x2,y2,z2);
	return leth;
}
double leastsquare(double *x,double *y,double *p,double n)//最小二乘法拟合自动去掉y<=0的点，返回r值
{
	double r,xy,xb,yb,xx,t,yy;
	t=n;
	r=xy=xb=yb=xx=yy=0;
	for(int i=0;i<n;i++)
	{
		if(y[i]<=0)
		{
			t--;
			continue;
		}
		xb+=x[i];
		yb+=y[i];
		xy+=x[i]*y[i];
		xx+=x[i]*x[i];
		yy+=y[i]*y[i];
	}
	xb=xb/t;
	yb=yb/t;
	xy=xy/t;
	xx=xx/t;
	yy=yy/t;
	p[0]=(xy-xb*yb)/(xx-xb*xb);
	p[1]=yb-p[0]*xb;
	r=(xy-xb*yb)/sqrt((xx-xb*xb)*(yy-yb*yb));
	//for(int i=0;i<n;i++)
	//{
	//sum1
	return r;
}
//给出x-z和y-z拟合结果p1和p2，组合出空间的直线方程，并求出给定点（x,y,z）到直线的距离
double dline(double *p1,double *p2,double x,double y)
{
	double d;
	//d=fabs(p2[0]*y-p1[0]*x+p2[1]-p1[1])/sqrt(p1[0]*p1[0]+p2[0]*p2[0]);
	d=fabs(p1[0]*x-p2[0]*y+p1[1]-p2[1])/sqrt(p1[0]*p1[0]+p2[0]*p2[0]);
	
	return d;
}
double dline1(double *p,double x,double y,double z)
{
	double a1,a2,a3,t,l1,l2,d;
	t=1;
	a1=p[0]+p[1]*t-x;
	a2=p[2]+p[3]*t-y;
	a3=t-z;
	l1=(p[1]*a1+p[3]*a2+a3)/sqrt(p[1]*p[1]+p[3]*p[3]+1);
	l2=sqrt(a1*a1+a2*a2+a3*a3);
	d=sqrt(l2*l2-l1*l1);
	return d;
}
//重建直线的轨迹，输入数据为对应的(x,a),(y1,a1),(y2,a2);K为第K个事件，有N[K]个数据点
void reconstruct(int *xl,int *y1,int *y2,int *al,int *a1,int *a2,int K)
{
	N[K]=0;//初始设定N[K]为0,以后每对上一个点就++
	double error=6;
	double yy[80],aa[80],x[48],a[48];
	for(int i=0;i<80;i++) //首先合并y1和y2两个数组
	{
		if(i<48)
		{
			x[i]=xl[i];
			a[i]=al[i];
		}
		if(i<40)
		{
			yy[i]=y1[i];
			aa[i]=a1[i];
		}
		if(i>=40)
		{
			yy[i]=y2[i-40];
			aa[i]=a2[i-40];
		}
	}
	double tx[48],ty[80];
	for(int i=0;i<80;i++)
	{
		if(i<48)
		{
			tx[i]=(i+1)*2.08333;//换算成坐标
			x[i]=x[i]/10;//近似漂移速度相等
		}
		if(i<40)
		{
			ty[i]=(i+1)*2.5;
			yy[i]=(yy[i]+resmean[i])*velocity[i]/100;
			//yy[i]=yy[i]/10;
		}
		if(i>=40)
		{
			ty[i]=(i-39)*2.5;
			yy[i]=(yy[i]+resmean[i])*velocity[i]/100;
			//yy[i]=yy[i]/10;
		}
	}

	double par1[2],par2[2];
	leastsquare(tx,x,par1,48);
	leastsquare(ty,yy,par2,80);
	//cout<<"rx="<<leastsquare(tx,x,par1,48)<<"  k="<<par1[0]<<"  b="<<par1[1]<<endl;
	//cout<<"ry="<<leastsquare(ty,yy,par2,80)<<"  k="<<par2[0]<<"  b="<<par2[1]<<endl;
	double pl[4]={-par1[1]/par1[0],1/par1[0],-par2[1]/par2[0],1/par2[0]};

//x->y
	for(int i=0;i<48;i++)
	{
		if(x[i]<=0)  //把没有意义的点去掉
		continue;
		double smin1,smin2,stan;
		int t1,t2;
		smin1=smin2=9999;//一个在左边一个在右边
		t1=t2=-1;
		for(int j=0;j<80;j++)
		{
			if(yy[j]<=0)//把没有意义的点去掉
			continue;
			stan=dline1(pl,tx[i],ty[j],x[i]);
			//stan=fabs(yy[j]-x[i])+dline(pl,tx[i],ty[j],x[i]);
			//stan=fabs(yy[j]-x[i]);
			if(yy[j]==x[i])
			{
				t1=t2=j;
				smin1=smin2=0;
			}
			if(yy[j]<x[i])
			{
				if(stan<smin1)
				{
					smin1=stan;
					t1=j;
				}
			}
			else
			{
				if(stan<smin2)
				{
					smin2=stan;
					t2=j;
				}
			}
		}
		if(t1*t2<0)
		{
			//cout<<"没有找到对应的两点"<<endl;
			continue;
		}
		if(fabs(smin1)>error)//如果离它最近的点都距离很大，直接舍去这个点
		{
			//cout<<"xsmin距离太大"<<endl;
			continue;
		}
		if(fabs(yy[t1]-yy[t2])>error)
		{
			//cout<<"y1和y2距离过大"<<endl;
			continue;
		}
		else
		{	
			if(yy[t2]==yy[t1])
			{
				Y[K][N[K]]=ty[t1];
			}
			else
				Y[K][N[K]]=(x[i]-yy[t1])*(ty[t2]-ty[t1])/(yy[t2]-yy[t1])+ty[t1];
			//cout<<"t1="<<t1<<"  t2="<<t2<<endl;
			X[K][N[K]]=tx[i];
			Z[K][N[K]]=x[i];
			A[K][N[K]]=a[i];
			N[K]++;
			//cout<<"yy2="<<yy[t2]<<"  yy1="<<yy[t1]<<endl;
		}
	}

//y->x
	for(int i=0;i<80;i++)
	{
		if(yy[i]<=0)
		continue;
		double smin1,smin2,stan;
		int t1,t2;
		smin1=smin2=9999;
		t1=t2=-1;
		for(int j=0;j<48;j++)
		{
			if(x[j]<=0)//把没有意义的点去掉
			continue;
			stan=dline1(pl,tx[j],ty[i],yy[i]);
			//stan=fabs(x[j]-yy[i])+1.7*dline(par1,par2,tx[j],ty[i]);
			//stan=fabs(x[j]-yy[i]);
			if(x[j]==yy[i])
			{
				t1=t2=j;
				smin1=smin2=0;
			}
			if(x[j]<yy[i])
			{
				if(stan<smin1)
				{
					smin1=stan;
					t1=j;
				}
			}
			else
			{
				if(stan<smin2)
				{
					smin2=stan;
					t2=j;
				}
			}
		}
		//cout<<"tx1="<<tx[t1]<<" x1="<<x[t1]<<" tx2="<<tx[t2]<<" x2="<<x[t2]<<endl;
		//cout<<"y[i]="<<yy[i]<<" ty="<<ty[i]<<" smin1="<<smin1<<endl;
		if(t1*t2<0)
		{
			//cout<<"没有找到对应的两点"<<endl;
			continue;
		}
		if(fabs(smin1)>error)
		{
			//cout<<"ysmin太大"<<endl;
			continue;
		}
		if(fabs(x[t2]-x[t1])>error)
		{
			//cout<<"x1和x2距离过大"<<endl;
			continue;
		}
		else
		{
			if(x[t2]==x[t1])
				X[K][N[K]]=tx[t1];
			else
				X[K][N[K]]=(yy[i]-x[t1])*(tx[t2]-tx[t1])/(x[t2]-x[t1])+tx[t1];
			//cout<<"t1="<<t1<<"  t2="<<t2<<endl;
			Y[K][N[K]]=ty[i];
			Z[K][N[K]]=yy[i];
			A[K][N[K]]=aa[i];
			N[K]++;
			//cout<<"x1="<<x[t1]<<"  x2="<<x[t2]<<endl;
		}
	}
}
void reconstruct2(int *xl,int *y1,int *y2,int *al,int *a1,int *a2,int K)
{
	N[K]=0;//初始设定N[K]为0,以后每对上一个点就++
	double error=10;
	double yy[80],aa[80],x[48],a[48];
	for(int i=0;i<80;i++) //首先合并y1和y2两个数组
	{
		if(i<48)
		{
			a[i]=al[i];
			x[i]=xl[i];
			//if(a[i]>=0)
			//x[i]+=mc[K];
			
		}
		if(i<40)
		{
			yy[i]=y1[i];
			aa[i]=a1[i];
		}
		if(i>=40)
		{
			yy[i]=y2[i-40];
			aa[i]=a2[i-40];
		}
	}
	double tx[48],ty[80];
	for(int i=0;i<80;i++)
	{
		if(i<48)
		{
			tx[i]=(i+1)*2.08333;//换算成坐标
			x[i]=x[i]/10;//近似漂移速度相等
		}
		if(i<40)
		{
			ty[i]=(i+1)*2.5;
			yy[i]=yy[i]/10;
			//yy[i]=(yy[i])*velocity[i]/100;
		}
		if(i>=40)
		{
			ty[i]=(i-39)*2.5;
			yy[i]=yy[i]/10;
			//yy[i]=(yy[i])*velocity[i]/100;
		}
	}
	double yy1[40],yy2[40],ty1[40],ty2[40];
	int p,q;
	p=q=0;
	for(int i=0;i<80;i++)
	{
		if(ty[i]<=50)
		{
			ty1[p]=ty[i];
			yy1[p]=yy[i];
			p++;
		}
		else
		{
			ty2[q]=ty[i];
			yy2[q]=yy[i];
			q++;
		}
	}
	double par1[2],par2[2],par3[2];
	leastsquare(tx,x,par1,48);
	leastsquare(ty1,yy1,par2,40);
	leastsquare(ty2,yy2,par3,40);
	//cout<<"rx="<<leastsquare(tx,x,par1,48)<<"  k="<<par1[0]<<"  b="<<par1[1]<<endl;
	//cout<<"ry1="<<leastsquare(ty1,yy1,par2,40)<<"  k="<<par2[0]<<"  b="<<par2[1]<<endl;
	//cout<<"ry2="<<leastsquare(ty2,yy2,par3,40)<<"  k="<<par3[0]<<"  b="<<par3[1]<<endl;
	double pl[4]={-par1[1]/par1[0],1/par1[0],-par2[1]/par2[0],1/par2[0]};
	//Drawline(pl);
//x->y
	for(int i=0;i<48;i++)
	{
		if(x[i]<=0)  //把没有意义的点去掉
		continue;
		if(x[i]<=((yy[19]>yy[59])?yy[19]:yy[59]))
		{
			if((x[i]/par2[0]-par2[1]/par2[0])>100||(x[i]/par2[0]-par2[1]/par2[0])<0)
			{
				//cout<<"x->y越界"<<endl;
				continue;
			}
			X[K][N[K]]=tx[i];
			Y[K][N[K]]=x[i]/par2[0]-par2[1]/par2[0];
			Z[K][N[K]]=x[i];
			A[K][N[K]]=a[i];
			//cout<<"下半段 "<<X[K][N[K]]<<" "<<Y[K][N[K]]<<" "<<Z[K][N[K]]<<endl;
			N[K]++;
			
		}
		else
		{
			if((x[i]/par3[0]-par3[1]/par3[0])>100||(x[i]/par3[0]-par3[1]/par3[0])<0)
			{
				//cout<<"x->y越界"<<endl;
				continue;
			}
			X[K][N[K]]=tx[i];
			Y[K][N[K]]=x[i]/par3[0]-par3[1]/par3[0];
			Z[K][N[K]]=x[i];
			A[K][N[K]]=a[i];
			//cout<<"上半段 "<<X[K][N[K]]<<" "<<Y[K][N[K]]<<" "<<Z[K][N[K]]<<endl;
			N[K]++;
			
		}
	}
//y->x
	for(int i=0;i<80;i++)
	{
		if(yy[i]<0)
		{
			continue;
		}
		if((yy[i]/par1[0]-par1[1]/par1[0])>100||(yy[i]/par1[0]-par1[1]/par1[0])<0)
		{
			//cout<<"y->x越界"<<endl;
			continue;
		}
		Y[K][N[K]]=ty[i];
		X[K][N[K]]=yy[i]/par1[0]-par1[1]/par1[0];
		Z[K][N[K]]=yy[i];
		A[K][N[K]]=aa[i];
		//cout<<"y 部分 "<<X[K][N[K]]<<" "<<Y[K][N[K]]<<" "<<Z[K][N[K]]<<endl;
		N[K]++;
		
	}
}
				



int mainline(int nstart)
{
	
	for(int s=0;s<event;s++)
	{N[s]=0;}

	int index_branch=-1;  //第几个筛选出来的事件
	int switch_leaf=0;  //事件中的第几类数据
	ifstream infile;
	infile.open("/home/jeff/TPC/data_331.txt"); //改成你的路径
	//string strTmp;
	long iTmp; //读到的数
	int m;
	for(int i=0; i<257*event;++i)
	{
		infile>>iTmp;
		if(iTmp<0) { iTmp=-1; }
		int jaziel = i%257; //每个事件包含257个数，jaziel记录该数在事件中的序号
		if(jaziel == 0)    //第一个数为事件编号
		{       m=0;
			++index_branch;
			switch_leaf=0;
			number[index_branch]=iTmp;
		}
		else if(jaziel<=40)switch_leaf=1; //第2~41个数为y1_tdc
		else if(jaziel<=80)switch_leaf=2; //第42~81个数为y2_tdc
		else if(jaziel<=120)switch_leaf=3;//第82~121个数为y1_adc
		else if(jaziel<=160)switch_leaf=4;//第122~161个数为y2_adc
		else if(jaziel<=208)switch_leaf=5;//第162~209个数为x_tdc
		else switch_leaf=6;               //第210~257个数为x_adc
		switch(switch_leaf)
		{
			case 0:
        			break;
      			case 1:
        			y01[index_branch][jaziel-1]=iTmp;
        			break;
      			case 2:
        			//y2[index_branch][jaziel-41]=iTmp;
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
				/*if(iTmp<2000)
				{
					m++;
					for(int j=jaziel-209;j<47;j++)
					{
						x[index_branch][j]=x[index_branch][j+1];
					}
					break;
				}
        			a3[index_branch][jaziel-209-m]=iTmp;
        			break;*/
				if(iTmp<2000)
				{
					x[index_branch][jaziel-209]=-1;
					a3[index_branch][jaziel-209]=-1;
					break;
				}
				a3[index_branch][jaziel-209]=iTmp;
				if(iTmp>5000)
				{
					x[index_branch][jaziel-209]+=-36*iTmp/(10000+iTmp)+12;//7.26896+1.69586
				}
				break;
      			default:
        			cout<<jaziel<<"  error!!"<<endl;

		}
		


	}
	SetHstyle();
	TCanvas *c1=new TCanvas("c1","",600,600);
  	TCanvas *c2=new TCanvas("c2","",600,600);
  	TCanvas *c3=new TCanvas("c3","",600,600);
	TH1F* Ga=new TH1F("Ga","Angle",100,0,5);
	TH1F* Gd=new TH1F("Gd","Distance",100,0,3);
	
	TH1F* Gl=new TH1F("Gl","Length",100,0,3);
																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																	    		TGraph2D *g=new TGraph2D();
	TCanvas *c4=new TCanvas("c4","",600,600);
	TCanvas *c5=new TCanvas("c5","",600,600);
	TCanvas *c6=new TCanvas("c6","",600,600);	

	
  	c4->SetMargin(0.2,0.2,0.1,0.1);
  	c5->Divide(2,2);

  	TH2F *hadcx = new TH2F("hadcx","adc vs X ch",48,-0.5,47.5,200,50,46050);
	hadcx->SetStats(false);
	TH2F *hadcy1 = new TH2F("hadcy1","adc vs Y1 ch",40,-0.5,39.5,200,50,46050);
	hadcy1->SetStats(false);
	//TH2F *hadcy2 = new TH2F("hadcy2","adc vs Y2 ch",40,39.5,79.5,200,1,101);
	//hadcy2->SetStats(false);
	//TH2F *htdcx = new TH2F("htdcx","tdc vs X ch",48,-0.5,47.5,200,-000,2000);
	TH2F *htdcx = new TH2F("htdcx","tdc vs X ch",48,0,100,200,0.1,140);
	htdcx->SetStats(false);
   	htdcx->SetMarkerStyle(5);
   	htdcx->SetMarkerSize(0.7);
   	htdcx->SetMarkerColor(9);
	//TH2F *htdcy1 = new TH2F("htdcy1","tdc vs Y1 ch",40,-0.5,39.5,200,-000,2000);
	TH2F *htdcy1 = new TH2F("htdcy1","tdc vs Y1 ch",40,0,100,200,0.1,140);
	htdcy1->SetStats(false);
   	htdcy1->SetMarkerStyle(5);
   	htdcy1->SetMarkerSize(0.7);
   	htdcy1->SetMarkerColor(9);
	//TH2F *htdcy2 = new TH2F("htdcy2","tdc vs Y2 ch",40,39.5,79.5,200,-000,2000);
	//htdcy2->SetStats(false);
   	//htdcy2->SetMarkerStyle(5);
   	//htdcy2->SetMarkerSize(0.7);
   	//htdcy2->SetMarkerColor(9);

  	double Yc[event],Xc[event];
	double ev=0;
	for(int k=0;k<event; ++k)
	{
		double p,q;
		p=q=0;
		Yc[k]=Xc[k]=mc[k]=0;
		for(int i=0;i<40;i++)
		{
			if(y01[k][i]<=0)
			{
				continue;
			}
			Yc[k]+=y01[k][i];
			p++;
		}
		for(int i=0;i<40;i++)
		{
			if(y2[k][i]<=0)
			{
				continue;
			}
			Yc[k]+=y2[k][i];
			p++;
		}
                Yc[k]=Yc[k]/p;
		for(int j=0;j<48;j++)
		{
			if(x[k][j]<=0)
			{
				continue;
			}
			Xc[k]+=x[k][j];
			q++;
		}
		if(p==0||q==0)
		{
			ev++;
			continue;
		}
		Xc[k]=Xc[k]/q;
		mc[k]=Yc[k]-Xc[k];
		//cout<<mc[k]<<endl;
	}
/*	double y_x;
	y_x=0;
	for(int k=0;k<event;k++)
	{
		y_x+=mc[k];
	}
	cout<<"yb-xb="<<y_x/(event-ev)<<endl;
*/



  	for(int k=nstart; k<event; ++k)
	{
    		double X2[120],Y2[120],Z2[120],X3[120],Y3[120],Z3[120];

    		int p=0,q=0,n=0;

    		//reconstruct1(number[k],x[k],y01[k],a1[k],a3[k],k);

    		//reconstruct1(number[k],x[k],y2[k],a2[k],a3[k],k);

		reconstruct2(x[k],y01[k],y2[k],a3[k],a1[k],a2[k],k);

    		Sort(X[k],Y[k],Z[k],A[k]);     //把点按照y从大到小的顺序排列，即粒子前进的方向

    
    		findpoint(N[k],X[k],Y[k],Z[k],A[k],X2,Y2,Z2,X3,Y3,Z3,&p,&q); //按平面把轨迹截断，同时去掉adc<1000的，分开的两端存在X2,X3中

		//findpoint1(N[k],X[k],Y[k],Z[k],A[k],X2,Y2,Z2,X3,Y3,Z3,&p,&q);//把轨迹按隔一个点分开
    		if(k%1000==0)
		cout<<"constructing "<<k<<"........"<<endl;

    		//cout<<number[k]<<endl;

  		g->SetTitle(";x/mm;y/mm;z/mm");

    		//g->SetMaximum(1000);



    		g->SetPoint(0,0,0,0);

    		g->SetPoint(1,100,100,140);  

		for(int j=0;j<p;j++)     //分开后前半段（后半段为q）的点
		{
			//cout<<X2[j]<<"\t"<<Y2[j]<<"\t"<<Z2[j]<<"\t"<<endl;
			g->SetPoint(j+2,X2[j],Y2[j],Z2[j]);
			//g->SetPoint(j+2,X3[j],Y3[j],Z3[j]);
		} 
		for(int j=0;j<q;j++)     //分开后前半段（后半段为q）的点
		{
			//cout<<X2[j]<<"\t"<<Y2[j]<<"\t"<<Z2[j]<<"\t"<<endl;
			g->SetPoint(j+2+p,X3[j],Y3[j],Z3[j]);
			//g->SetPoint(j+2,X3[j],Y3[j],Z3[j]);
		}             

 		/*for(int m=0;m<N[k];m++)  //第k组所有的点

    		{
         
        		cout<<X[k][m]<<"\t"<<Y[k][m]<<"\t"<<Z[k][m]<<"\t"<<A[k][m]<<endl;

       			g->SetPoint(m+2,X[k][m],Y[k][m],Z[k][m]);

		

    		} */

               c4->cd();
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
                      
		if(line3Dfit(X2,Z2,Y2,p,par1,ps1)==1) //拟合左半部分散点
			continue;
		DrawFitline1 (par1);          //画出左半部分的拟合线
	
	
	
		if(line3Dfit(X3,Z3,Y3,q,par2,ps2)==1) //拟合右半部分散点     
			continue;                
		DrawFitline2 (par2);          //画出右半部分的拟合线
	
		d[k]=footpoint(par1,par2,M);      //求出拟合线公垂线的中点，存在M[3]中
  		Gd->Fill(d[k]);

		a[k]=angle(par1,par2);            //角度
				
		Ga->Fill(a[k]);

		leth[k]=length(par1,par2);
		Gl->Fill(leth[k]);
		cout<<"angle="<<a[k]<<endl;
		cout<<"交点坐标为("<<M[0]<<","<<M[1]<<","<<M[2]<<")"<<endl;
		cout<<"两直线之间距离为"<<d[k]<<endl;
		cout<<"两直线分界面的距离为"<<leth[k]<<endl;

	        TGraph2D *g1=new TGraph2D();
		g1->SetPoint(0,0,0,0);
	
    		g1->SetPoint(1,100,100,140);
		for(int m=0;m<N[k];m++)  //第k组所有的点

    		{
         
        		//cout<<X[k][m]<<"\t"<<Y[k][m]<<"\t"<<Z[k][m]<<"\t"<<A[k][m]<<endl;

       			g1->SetPoint(m+2,X[k][m],Y[k][m],Z[k][m]);

		

    		}

		//g->SetPoint(100,M[0],M[1],M[2]);
	
		c6->cd();	
		g1->Draw(); 

 		for(int l=0;l<40;l++)
		{
	
			hadcy1->Fill(l,a1[k][l]);
			//htdcy1->Fill(l,tdc_y[l]);
			hadcy1->Fill(l,a2[k][l]);
			htdcy1->Fill(float(l)*2.5,y01[k][l]/10);
			htdcy1->Fill(float(l)*2.5,y2[k][l+40]/10);
		}
		for(int j=0;j<48;j++)
		{
			hadcx->Fill(j,a3[k][j]);
			htdcx->Fill(float(j)*2.08333,x[k][j]/10);
			//cout<<"xtdc["<<j<<"]="<<x[k][j]<<"  xadc["<<j<<"]="<<a3[k][j]<<endl;

		}
		gStyle->SetPalette(100);
        	c5->cd(1);
		hadcx->Draw("colz");
        	c5->cd(2);
		hadcy1->Draw("colz");  
       		c5->cd(3);
		htdcx->Draw("colz");  
        	c5->cd(4);
		htdcy1->Draw("colz");                                                        
		               break;                      
    		//c1->cd();

                //g->Draw();                                                             

		
        }
	double v[100];
	double sum2,sum1;
	sum2=sum1=0;
	for(int i=1;i<=100;i++)
	{
		v[i]=Ga->GetBinContent(i);
		//cout<<v[i]<<endl;
		v[i]=v[i]/sin((i*0.05-0.025)*3.1416/180);
		sum2+=((i*0.05-0.025)*(i*0.05-0.025))*v[i];
		sum1+=v[i];
		Ga->SetBinContent(i,v[i]);
	}
	cout<<"Angle_mean_sqr_root="<<sqrt(sum2/sum1)<<endl;
	sum2=sum1=0;
	for(int i=1;i<=100;i++)
	{
		v[i]=Gl->GetBinContent(i);
		sum2+=((i*0.03-0.015)*(i*0.03-0.015))*v[i];
		sum1+=v[i];
	}
	cout<<"Length_mean_sqr_root="<<sqrt(sum2/sum1)<<endl;
		
	gStyle->SetOptStat("nemr"); // the default value
	c1->cd();
	Ga->Draw();
	c2->cd();
	Gd->Draw();
	c3->cd();
	Gl->Draw();
return 0;
}

