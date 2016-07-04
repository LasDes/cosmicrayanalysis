/*Code to plot traces of muons detected in a single module of GRAPES-3
**Sanket Deshpande**
*/

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
 
#include <TBox.h>
#include <TLine.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <dirent.h>
#include "GROOTHeaders.h"
#include "GSetLeafAddress.h"
#include "GDateTimeHandler.h"
#include "GControlFileHandler.h"
#include "GCorsikaROOTReader.h"

using namespace std;

ofstream  errlog;


int main ()
{
	
	//reading the data file
	long int highest=0;	
	struct DateTimeHandler_Struct  sDTH;
	GDateTimeHandler cDTH;

	GCorsikaROOTReader aa("../DAT000020.root");
	if(aa.IsInputOK()) {exit(0);};
	
	int nentries = aa.GetEntries();
		
	for(long int i=0; i<nentries; i++) {
		aa.GetEntry(i);
		if(aa.np>highest)
			highest=aa.np;
	}
	//cout<<highest<<endl;

	float xmin[232];
	float xmax[232];
	float ymin[232];
	float ymax[232];
	float zmin[232];
	float zmax[232];
	float xmaxmod=0.0, xminmod=0.0,ymaxmod=0.0, yminmod=0.0, centrex, centrey;
	int n=0;

	int station, module, layer[232], counter[232];	

	//reading the prc file	
	ifstream fileprc;
	fileprc.open("prc_coordinates.txt");
	while(fileprc>>station)
	{
		fileprc>>module;
		if(module>0)
			break;
		fileprc>>layer[n];
		fileprc>>counter[n];
		fileprc>>xmin[n];
		fileprc>>xmax[n];
		fileprc>>ymin[n];
		fileprc>>ymax[n];
		fileprc>>zmin[n];
		fileprc>>zmax[n];
		n++;
	}

	//finding the corner points of the module
	for(int i=0; i<=n; i++)
	{
		if(xmax[i]>xmaxmod)
			xmaxmod=xmax[i];

		if(xmin[i]<xminmod)
			xminmod=xmin[i];
	
		if(ymax[i]>ymaxmod)
			ymaxmod=ymax[i];
	
		if(ymin[i]<yminmod)
			yminmod=ymin[i];
	}


	//coordinates of the centre of the module
	centrex=(xmaxmod+xminmod)/2;
	centrey=(ymaxmod+yminmod)/2;

	//shifting the entire module to the centre of the shower core	
	for (int i=0; i<=n; i++)
	{
		xmax[i]=xmax[i]-centrex;
		xmin[i]=xmin[i]-centrex;
		ymax[i]=ymax[i]-centrey;
		ymin[i]=ymin[i]-centrey;
	}
	float a, b, c, d;
	TCanvas *c1[nentries];
	
	TBox box1[2][58];
	TBox box2[2][58];
	TFile *f = new TFile("Muon_Trace.root","RECREATE");

	for (int shower=0; shower<nentries; shower++)
	{
		aa.GetEntry(shower);


		int zdiffone, zdifftwo, flag;
		float mom, anglez, anglex, rone, rtwo;
		int layer3counter[4], layer2counter[4], layer1counter[4], layer0counter[4], detection[4][58];
		int layer3det[4], layer2det[4], layer1det[4], layer0det[4];
		int ctt=0, ct=0, co=0, cz=0;
		for(int i=0; i<=3; i++)
		{
			for(int j=0; j<58; j++)
				detection[i][j]=0;
		}

		for(int i=0; i<4; i++)
		{
			layer3counter[i]=60;
			layer2counter[i]=60;
			layer1counter[i]=60;
			layer0counter[i]=60;
			layer3det[i]=60;
			layer2det[i]=60;
			layer1det[i]=60;
			layer0det[i]=60;
		}

		int trace1[aa.np], trace2[aa.np];
		int det1=0, det2=0;
		int particleno1[aa.np], particleno2[aa.np];

		for(int i=0; i<aa.np; i++)
		{
			trace1[i]=0;
			trace2[i]=0;
		}

		float Xeffone, Xefftwo, Yeffone, Yefftwo;

		//if (ctt>0 || ct>0 || co>0 || cz>0)
			cout<<"\nSHOWER #: "<<shower<<endl;
		for(int j=0; j<aa.np; j++)
		{ 
			if (aa.PId[j]==5 || aa.PId[j]==6)
			{
				mom=sqrt(aa.Px[j]*aa.Px[j]+aa.Py[j]*aa.Py[j]+aa.Pz[j]*aa.Pz[j]);//total momentum
				anglez=acos(aa.Pz[j]/mom); //angle made by the net momentum with the z axis
				anglex=acos(aa.Px[j]/(mom*sin(anglez))); //angle made by the projection of net momentum on the x-y plane with the x axis

				for(int b=0; b<=231; b++)
				{
	
					switch (layer[b])
					{
						//for layer 3						
						case 3:
				 		flag=0;
						zdiffone=0; //differnce of this layer's zmax and upper layer's zmin
						zdifftwo=10; //difference of the layer's zmin and upper layer's zmin
						rone=zdiffone*tan(anglez);
						rtwo=zdifftwo*tan(anglez);
						Xeffone= aa.X[j] + rone*cos(anglex);
						Xefftwo= aa.X[j] + rtwo*cos(anglex);
						Yeffone= aa.Y[j] + rone*sin(anglex);
						Yefftwo= aa.Y[j] + rtwo*sin(anglex);
						//checking the track of the particle, both at the top of the counter and at the bottom
						//*one[][] represent the top whereas *two represent the bottom of the counter 
						if(Xeffone>=xmin[b] && Xeffone<=xmax[b] && Yeffone>=ymin[b] && Yeffone<=ymax[b])
						{				
							layer3counter[ctt]=counter[b]; //recording the counter number
							cout<<"Layer 3, counter= "<<layer3counter[ctt]<<" particle # "<<j<<endl;
							detection[3][counter[b]]++; //marking the particular counter of the layer 3 of having detected a muon
							trace1[j]++; 
							ctt++;
							flag++;
						}
						if(Xeffone>=xmin[b] && Xeffone<=xmax[b] && Yeffone>=ymin[b] && Yeffone<=ymax[b] && flag==0)
						{				
							layer3counter[ctt]=counter[b];
							cout<<"Layer 3, counter= "<<layer3counter[ctt]<<" particle # "<<j<<endl;
							detection[3][counter[b]]++; //marking the particular counter of the layer 3 of having detected a muon
							trace1[j]++;
							ctt++;
						}
						break;

						//for layer 2
						case 2:
						flag=0;
						zdiffone=10+17; //differnce of this layer's zmax and upper layer's zmin
						zdifftwo=10+17+10; //difference of the layer's zmin and upper layer's zmin
						rone=zdiffone*tan(anglez);
						rtwo=zdifftwo*tan(anglez);
						Xeffone= aa.X[j] + rone*cos(anglex); //Xeff at the top of the counter
						Xefftwo= aa.X[j] + rtwo*cos(anglex); //Xeff at the bottom of the counter
						Yeffone= aa.Y[j] + rone*sin(anglex); //
						Yefftwo= aa.Y[j] + rtwo*sin(anglex);
						if(Xeffone>=xmin[b] && Xeffone<=xmax[b] && Yeffone>=ymin[b] && Yeffone<=ymax[b])
						{				
							layer2counter[ct]=counter[b];
							cout<<"Layer 2, counter= "<<layer2counter[ct]<<" particle # "<<j<<endl;
							detection[2][counter[b]]++; //marking a particular counter of the layer 2 of having detected a muon
							trace2[j]++;
							ct++;
							flag++;
						}
						if(Xeffone>=xmin[b] && Xeffone<=xmax[b] && Yeffone>=ymin[b] && Yeffone<=ymax[b] && flag==0)
						{				
							layer2counter[ct]=counter[b];
							cout<<"Layer 2, counter= "<<layer2counter[ct]<<" particle # "<<j<<endl;
							detection[2][counter[b]]++; //marking a particular counter of the layer 2 of having detected a muon
							trace2[j]++;
							ct++;
						}
						break;

						//for layer 1
						case 1:
						flag=0;
						zdiffone=10+17+10+16;
						zdifftwo=10+17+10+16+10;
						rone=zdiffone*tan(anglez);
						rtwo=zdifftwo*tan(anglez);
						Xeffone= aa.X[j] + rone*cos(anglex);
						Xefftwo= aa.X[j] + rtwo*cos(anglex);
						Yeffone= aa.Y[j] + rone*sin(anglex);
						Yefftwo= aa.Y[j] + rtwo*sin(anglex);
						if(Xeffone>=xmin[b] && Xeffone<=xmax[b] && Yeffone>=ymin[b] && Yeffone<=ymax[b])
						{				
							layer1counter[co]=counter[b];
							cout<<"Layer 1, counter= "<<layer1counter[co]<<" particle # "<<j<<endl;
							detection[1][counter[b]]++;
							trace1[j]++;
							co++;
							flag++;
						}
						if(Xeffone>=xmin[b] && Xeffone<=xmax[b] && Yeffone>=ymin[b] && Yeffone<=ymax[b] && flag==0)
						{	layer1counter[co]=counter[b];
							cout<<"Layer 1, counter= "<<layer1counter[co]<<" particle # "<<j<<endl;
							detection[1][counter[b]]++;
							trace1[j]++;
							co++;
						}
						break;

						//for layer 0
						case 0:
						flag=0;
						zdiffone=10+17+10+16+10+16;
						zdifftwo=10+17+10+16+10+16+10;
						rone=zdiffone*tan(anglez);
						rtwo=zdifftwo*tan(anglez);
						Xeffone= aa.X[j] + rone*cos(anglex);
						Xefftwo= aa.X[j] + rtwo*cos(anglex);
						Yeffone= aa.Y[j] + rone*sin(anglex);
						Yefftwo= aa.Y[j] + rtwo*sin(anglex);
						if(Xeffone>=xmin[b] && Xeffone<=xmax[b] && Yeffone>=ymin[b] && Yeffone<=ymax[b])
						{				
							layer0counter[cz]=counter[b];
							cout<<"Layer 0, counter= "<<layer0counter[cz]<<" particle # "<<j<<endl;
							detection[0][counter[b]]++;
							trace2[j]++;
							cz++;
							flag++;
						}
						if(Xeffone>=xmin[b] && Xeffone<=xmax[b] && Yeffone>=ymin[b] && Yeffone<=ymax[b] && flag==0)
						{				
							layer0counter[cz]=counter[b];
							cout<<"Layer 0, counter= "<<layer0counter[cz]<<" particle # "<<j<<endl;
							detection[0][counter[b]]++;
							trace2[j]++;							
							cz++;
						}
		
						break;
					}
		
				}
				//only those particles which are detected by both layer 1 and 3 can make trace1=2, hence, the following conditional statement eliminates those particles which do not allow us to trace their paths
				//the layer*det[] arrays store the counter numbers of respective layers which have been struck by the same muon, in a sequential fashion
				if (trace1[j]>=2)
				{
					
					layer3det[det1]=layer3counter[ctt-1];
					layer1det[det1]=layer1counter[co-1];
					particleno1[det1]=j;
					det1++;
					
				}
				
				
				//only those particles which are detected by both layer 0 and 2 can make trace2=2, hence, the following conditional statement eliminates those particles which do not allow us to trace their paths	
				if (trace2[j]>=2)
				{
					layer2det[det2]=layer2counter[ct-1];
					layer0det[det2]=layer0counter[cz-1];
					particleno2[det2]=j;
					det2++;
				}

					
		
				if (trace1[j]>=2 || trace2[j]>=2)
				{
					cout<<"\nMuon path tracing \n";
					if(det1!=0)
					{
					for(int i=0; i<det1; i++)
						cout<<"\tParticle no "<<particleno1[i]<<" struck layer1 counter # "<<layer1det[i]<<" and layer3 counter # "<<layer3det[i]<<endl;
					}
					else
						cout<<"\nCannot trace for layer 1 and 3"<<endl;

					if(det2!=0)
					{
					for(int i=0; i<det2; i++)
						cout<<"\tParticle no "<<particleno2[i]<<" struck layer0 counter # "<<layer0det[i]<<" and layer2 counter # "<<layer2det[i]<<endl;
					}
					else
						cout<<"\nCannot trace for layer 0 and 2"<<endl;
				}
			}
		}
		
		
		if (ctt>0)
		{
			cout<<"\n For layer 3, following counters recorded muons: "<<endl;
			for(int i=0; i<ctt; i++)
			{
				if(layer3counter[i]!=60)
					cout<<layer3counter[i]<<endl;	
			}
		}

		if (ct>0)
		{
			cout<<"\n For layer 2, following counters recorded muons: "<<endl;
			for(int i=0; i<ct; i++)
			{
				if (layer2counter[i]!=60)
					cout<<layer2counter[i]<<endl;
			}
		}

		if (co>0)
		{
			cout<<"\n For layer 1, following counters recorded muons: "<<endl;
			for(int i=0; i<co; i++)
			{
				if (layer1counter[i]!=60)		
					cout<<layer1counter[i]<<endl;
			}
		}

		if (cz>0)
		{
			cout<<"\n For layer 0, following counters recorded muons: "<<endl;
			for(int i=0; i<cz; i++)	
			{
				if (layer0counter[i]!=60)
					cout<<layer0counter[i]<<endl;
			}
		}
		
		c1[shower] = new TCanvas(Form("c%d",shower));
		for(int i=0; i<4; i=i+2)
		{	
			b=0.2 + i/10.0;
			d=0.225 + i/10.0;
		
			for(int k=1; k<59; k++)
			{	a=k/59.0;
				c=(k+0.8)/59.0;
			
				if(detection[i][k-1]>0)
				{
					if (i==2)
					{
						box1[i-2][k-1].SetFillColor(2);
						box1[i-2][k-1].DrawBox( a,b,c,d);
						box2[i-2][k-1].SetFillStyle(0);
	   					box2[i-2][k-1].SetLineColor(1);
	   					box2[i-2][k-1].SetLineWidth(1);
						box2[i-2][k-1].DrawBox( a,b,c,d);
					}
					else
					{
						box1[i+1][k-1].SetFillColor(2);
						box1[i+1][k-1].DrawBox( a,b,c,d);
						box2[i+1][k-1].SetFillStyle(0);
	   					box2[i+1][k-1].SetLineColor(1);
	   					box2[i+1][k-1].SetLineWidth(1);
						box2[i+1][k-1].DrawBox( a,b,c,d);
					}
				}
				else
				{	
					if (i==2)
					{
						box2[i-1][k-1].SetFillStyle(0);
		   				box2[i-1][k-1].SetLineColor(1);
	   					box2[i-1][k-1].SetLineWidth(1);
						box2[i-1][k-1].DrawBox( a,b,c,d);
					}
					else
					{
						box2[i][k-1].SetFillStyle(0);
		   				box2[i][k-1].SetLineColor(1);
	   					box2[i][k-1].SetLineWidth(1);
						box2[i][k-1].DrawBox( a,b,c,d);
					}
				}
			}
		}

		for(int i=1; i<4; i=i+2)
		{	b=0.6 + i/10.0;
			d=0.625 + i/10.0;
		
			for(int k=1; k<59; k++)
			{	a=k/59.0;
				c=(k+0.8)/59.0;
			
				if(detection[i][k-1]>0)
				{
					if(i==1)
					{
						box1[i-1][k-1].SetFillColor(2);
						box1[i-1][k-1].DrawBox( a,b,c,d);
						box2[i-1][k-1].SetFillStyle(0);
	   					box2[i-1][k-1].SetLineColor(1);
	   					box2[i-1][k-1].SetLineWidth(1);
						box2[i-1][k-1].DrawBox( a,b,c,d);
					}
					else if(i==3)
					{
						box1[i-2][k-1].SetFillColor(2);
						box1[i-2][k-1].DrawBox( a,b,c,d);
						box2[i-2][k-1].SetFillStyle(0);
	   					box2[i-2][k-1].SetLineColor(1);
	   					box2[i-2][k-1].SetLineWidth(1);
						box2[i-2][k-1].DrawBox( a,b,c,d);
				
					}
				}
				else
				{	
					if(i==1)
					{
						box2[i-1][k-1].SetFillStyle(0);
		   				box2[i-1][k-1].SetLineColor(1);
	   					box2[i-1][k-1].SetLineWidth(1);
						box2[i-1][k-1].DrawBox( a,b,c,d);
					}
					else if(i==3)
					{
						box2[i-2][k-1].SetFillStyle(0);
		   				box2[i-2][k-1].SetLineColor(1);
	   					box2[i-2][k-1].SetLineWidth(1);
						box2[i-2][k-1].DrawBox( a,b,c,d);
					}
				
				}
			
			
			}
		}
		

		//displaying the Tracing details of the muon which has been recorded by both the sets of layers 
	
		/*
		TLine *line1 = new TLine();

		//Calculating angle using known data and then drawing reference traces
		cout<<"\nAngles made by Traces:"<<endl;
		
		float actualangle1[58], actualangle2[58];
		cout<<"For Layers 3 and 1\n"<<endl;
		double x1, y1, x2, y2;
		for (int i=0; i<det1; i++)
		{
				x1=(layer1det[i]+0.8)/59.0;
				x2=(layer3det[i]+0.8)/59.0;
				y1=0.6 + (1/10.0);
				y2=0.625 + (3/10.0);
				actualangle1[i]=atan((y2-y1)/(x2-x1));
				cout<<"For trace between counter# "<<layer3det[i]<<" and counter# "<<layer1det[i]<<" , the angle is "<<actualangle1[i]<<endl;
				
				line1->SetLineWidth(2);
				line1->SetLineColor(5);
				line1->Draw();
				line1->DrawLine(x1,y1,x2,y2);
			
			
		}		

		//drawing the traces
		cout<<"For Layers 2 and 0\n"<<endl;
		for (int i=0; i<det2; i++)
		{
				//cout<<"i= "<<i<<endl;				
				x1=(layer0det[i]+0.8)/59.0;
				x2=(layer2det[i]+0.8)/59.0;
				y1=0.2 + (0/10.0);
				y2=0.225 + (2/10.0);
				actualangle2[i]=atan((y2-y1)/(x2-x1));
				cout<<"For trace between counter# "<<layer2det[i]<<" and counter# "<<layer0det[i]<<" , the angle is "<<actualangle2[i]<<endl;
				
				line1->SetLineWidth(2);
				line1->SetLineColor(5);
				line1->Draw();
				line1->DrawLine(x1,y1,x2,y2);
			
		}
		
		*/		

		//trying to plot the trace using the counter hits
		
		int probabletrace[4][60][60]; //stores the counter number of each counter from the lower layer which is 'hit'; first dimension is for layer number, second for the counter # for upper layer, third for counter # of lower layer which recieves a 'hit'
		float xb[4][60], yb[4][60]; //to store the coordinates of each of the counters from the upper layer which are 'hit'
		float xa[4][60][60], ya[4][60][60]; //arrays to store the coordinates of each of the counters from the lower layer which are 'hit'
		float surexa[4][60], sureya[4][60];
		float anglearray[4][60][60]; //array to store the angle made by each 
		int serno[2][60];
		int surecounter[4][60];
		int mostweight[2][60];
		int weight[2][60][60];
		int hitsonupper[4][60], hitsonlower[4][60];

		//for layers 0 and 2
		for(int i=0; i<60; i++)
		{
			mostweight[0][i]=0;
			hitsonupper[2][i]=0;
			hitsonlower[0][i]=0;
			for(int j=0; j<60; j++)
				weight[0][i][j]=0;
		}
		for (int i=0; i<2; i++)
		{
			for(int j=0; j<60; j++)
				serno[i][j]=0;
		}
		
		float sureangle[4][60];
		
		
		//for layer # 0 & 2
		for(int k=0; k<58; k++)
		{	
			if(detection[2][k]>0)
			{
				
				xb[2][k]=(k+0.8)/59.0;
				yb[2][k]=0.225 + (2/10.0);
				//cout<<"\n"<<k<<"\n"<<endl;
				for (int i=0; i<58; i++)
				{
					if (detection[0][i]>0)
						{
							probabletrace[0][k][serno[0][k]]=i;
							//cout<<probabletrace[0][k][serno[0][k]]<<endl;
							serno[0][k]++;
							hitsonlower[0][k]++;
						}
				}
				for (int i=0; i<serno[0][k]; i++)
				{
					xa[0][k][i]=(probabletrace[0][k][i]+0.8)/59.0;
					ya[0][k][i]=0.2 + (0/10.0);
					anglearray[0][k][i]= atan((yb[2][k]-ya[0][k][i])/(xb[2][k]-xa[0][k][i]));
					//cout<<anglearray[0][k][i]<<endl;
				}
				hitsonupper[2][k]++;		
			}
		}
		
		for(int k=0; k<58; k++)
		{
			if(detection[2][k]>0)
			{	
				//cout<<"\n"<<k<<"\n"<<endl;	
				for (int i=0; i<=serno[0][k]; i++)
				{
					for(int j=k+1; j<58; j++)
					{
						for (int h=0; h<serno[0][j];h++)
						{
							if(anglearray[0][j][h]<=anglearray[0][k][i]+0.10472 && anglearray[0][j][h]>=anglearray[0][k][i]-0.10472)
	{
		weight[0][k][i]++;
		if(anglearray[0][j][h]<=anglearray[0][k][i]+0.0872665 && anglearray[0][j][h]>=anglearray[0][k][i]-0.0872665)
		{
			weight[0][k][i]++;
			if(anglearray[0][j][h]<=anglearray[0][k][i]+0.0698132 && anglearray[0][j][h]>=anglearray[0][k][i]-0.0698132)
			{
				weight[0][k][i]++;
				if(anglearray[0][j][h]<=anglearray[0][k][i]+0.0523599 && anglearray[0][j][h]>=anglearray[0][k][i]-0.0523599)
				{
					weight[0][k][i]++;
					if(anglearray[0][j][h]<=anglearray[0][k][i]+0.0349066 && anglearray[0][j][h]>=anglearray[0][k][i]-0.0349066)
					{
						weight[0][k][i]++;
						if(anglearray[0][j][h]<=anglearray[0][k][i]+0.0174533 && anglearray[0][j][h]>=anglearray[1][k][i]-0.0174533)
						{
							weight[0][k][i]++;
							if(anglearray[0][j][h]<=anglearray[0][k][i]+0.00872665 && anglearray[0][j][h]>=anglearray[0][k][i]-0.00872665)
							{
								weight[0][k][i]++;
							}
						}
					}
				}
			}
		}
	}
						}
						
					}

					for(int j=k-1; j>=0; j--)
					{
						for (int h=0; h<serno[0][j];h++)
						{
							if(anglearray[0][j][h]<=anglearray[0][k][i]+0.10472 && anglearray[0][j][h]>=anglearray[0][k][i]-0.10472)
	{
		weight[0][k][i]++;
		if(anglearray[0][j][h]<=anglearray[0][k][i]+0.0872665 && anglearray[0][j][h]>=anglearray[0][k][i]-0.0872665)
		{
			weight[0][k][i]++;
			if(anglearray[0][j][h]<=anglearray[0][k][i]+0.0698132 && anglearray[0][j][h]>=anglearray[0][k][i]-0.0698132)
			{
				weight[0][k][i]++;
				if(anglearray[0][j][h]<=anglearray[0][k][i]+0.0523599 && anglearray[0][j][h]>=anglearray[0][k][i]-0.0523599)
				{
					weight[0][k][i]++;
					if(anglearray[0][j][h]<=anglearray[0][k][i]+0.0349066 && anglearray[0][j][h]>=anglearray[0][k][i]-0.0349066)
					{
						weight[0][k][i]++;
						if(anglearray[0][j][h]<=anglearray[0][k][i]+0.0174533 && anglearray[0][j][h]>=anglearray[1][k][i]-0.0174533)
						{
							weight[0][k][i]++;
							if(anglearray[0][j][h]<=anglearray[0][k][i]+0.00872665 && anglearray[0][j][h]>=anglearray[0][k][i]-0.00872665)
							{
								weight[0][k][i]++;
							}
						}
					}
				}
			}
		}
	}
						}
					}
					//cout<<weight[0][k][i]<<endl;				
				}
			}
		}
		

		cout<<"\n Reverse plotting the traces:"<<endl;
		
		cout<<" For layers 2 and 0"<<endl;
		for(int k=0; k<58; k++)
		{
		
			if(detection[2][k]>0)
			{
				for (int i=0; i<serno[0][k]; i++)
				{
					if(mostweight[0][k]<weight[0][k][i])
					{
						mostweight[0][k]=weight[0][k][i];
						sureangle[0][k]=anglearray[0][k][i];
						surecounter[0][k]=probabletrace[0][k][i];
						surexa[0][k]=xa[0][k][i];
						sureya[0][k]=ya[0][k][i];
					}
						
				}
				cout<<"mostweight for "<<k<<" is "<<mostweight[0][k]<<endl;
				if (mostweight[0][k]==0)
				{
					if(hitsonlower[0][k]==1)
					{	
						sureangle[0][k]=anglearray[0][k][0];
						surecounter[0][k]=probabletrace[0][k][0];
						surexa[0][k]=xa[0][k][0];
						sureya[0][k]=ya[0][k][0];
						goto display1;
					}
					else
					cout<<"There isnt a deciding value"<<endl;
					
				}
				else
				{
					display1:
					cout<<"Trace between: counter# "<<k<<" and counter # "<<surecounter[0][k]<<" has angle= "<<sureangle[0][k]<<endl;
					TLine *line2 = new TLine();
					line2->SetLineWidth(2);
					line2->SetLineColor(8);
					line2->Draw();
					line2->DrawLine(surexa[0][k],sureya[0][k],xb[2][k],yb[2][k]);
				}
				
			}
			

			
					
		}
		
		
		//for layers 1 and 3
		for(int i=0; i<60; i++)
		{
			mostweight[1][i]=0;
			hitsonupper[3][i]=0;
			hitsonlower[1][i]=0;
			for(int j=0; j<60; j++)
				weight[1][i][j]=0;
		}
		
		
	
		
		for(int k=0; k<58; k++)
		{	
			if(detection[3][k]>0)
			{
				
				xb[3][k]=(k+0.8)/59.0;
				yb[3][k]=0.625 + (3/10.0);
				//cout<<"\n"<<k<<"\n"<<endl;
				for (int i=0; i<58; i++)
				{
					if (detection[1][i]>0)
						{
							probabletrace[1][k][serno[1][k]]=i;
							//cout<<probabletrace[0][k][serno[0][k]]<<endl;
							serno[1][k]++;
							hitsonlower[1][k]++;
						}
				}
				
				for (int i=0; i<serno[1][k]; i++)
				{
					xa[1][k][i]=(probabletrace[1][k][i]+0.8)/59.0;
					ya[1][k][i]=0.6 + (1/10.0);
					anglearray[1][k][i]= atan((yb[3][k]-ya[1][k][i])/(xb[3][k]-xa[1][k][i]));
					//cout<<anglearray[0][k][i]<<endl;
				}
				hitsonupper[3][k]++;		
			}
		}
		
		for(int k=0; k<58; k++)
		{
			if(detection[3][k]>0)
			{	
				cout<<"\n"<<k<<"\n"<<endl;	
				for (int i=0; i<=serno[1][k]; i++)
				{
					for(int j=k+1; j<58; j++)
					{
						for (int h=0; h<serno[1][j];h++)
						{
							if(anglearray[1][j][h]<=anglearray[1][k][i]+0.10472 && anglearray[1][j][h]>=anglearray[1][k][i]-0.10472)
	{
		weight[1][k][i]++;
		if(anglearray[1][j][h]<=anglearray[1][k][i]+0.0872665 && anglearray[1][j][h]>=anglearray[1][k][i]-0.0872665)
		{
			weight[1][k][i]++;
			if(anglearray[1][j][h]<=anglearray[1][k][i]+0.0698132 && anglearray[1][j][h]>=anglearray[1][k][i]-0.0698132)
			{
				weight[1][k][i]++;
				if(anglearray[1][j][h]<=anglearray[1][k][i]+0.0523599 && anglearray[1][j][h]>=anglearray[1][k][i]-0.0523599)
				{
					weight[1][k][i]++;
					if(anglearray[1][j][h]<=anglearray[1][k][i]+0.0349066 && anglearray[1][j][h]>=anglearray[1][k][i]-0.0349066)
					{
						weight[1][k][i]++;
						if(anglearray[1][j][h]<=anglearray[1][k][i]+0.0174533 && anglearray[1][j][h]>=anglearray[1][k][i]-0.0174533)
						{
							weight[1][k][i]++;
							if(anglearray[1][j][h]<=anglearray[1][k][i]+0.00872665 && anglearray[1][j][h]>=anglearray[1][k][i]-0.00872665)
							{
								weight[1][k][i]++;
							}
						}
					}
				}
			}
		}					
	}
							
						}
						
					}

					for(int j=k-1; j>=0; j--)
					{
						for (int h=0; h<serno[1][j];h++)
						{
							if(anglearray[1][j][h]<=anglearray[1][k][i]+0.10472 && anglearray[1][j][h]>=anglearray[1][k][i]-0.10472)
	{
		weight[1][k][i]++;
		if(anglearray[1][j][h]<=anglearray[1][k][i]+0.0872665 && anglearray[1][j][h]>=anglearray[1][k][i]-0.0872665)
		{
			weight[1][k][i]++;
			if(anglearray[1][j][h]<=anglearray[1][k][i]+0.0698132 && anglearray[1][j][h]>=anglearray[1][k][i]-0.0698132)
			{
				weight[1][k][i]++;
				if(anglearray[1][j][h]<=anglearray[1][k][i]+0.0523599 && anglearray[1][j][h]>=anglearray[1][k][i]-0.0523599)
				{
					weight[1][k][i]++;
					if(anglearray[1][j][h]<=anglearray[1][k][i]+0.0349066 && anglearray[1][j][h]>=anglearray[1][k][i]-0.0349066)
					{
						weight[1][k][i]++;
						if(anglearray[1][j][h]<=anglearray[1][k][i]+0.0174533 && anglearray[1][j][h]>=anglearray[1][k][i]-0.0174533)
						{
							weight[1][k][i]++;
							if(anglearray[1][j][h]<=anglearray[1][k][i]+0.00872665 && anglearray[1][j][h]>=anglearray[1][k][i]-0.00872665)
							{
								weight[1][k][i]++;
							}
						}
					}
				}
			}
		}
	}
						}
					}
					//cout<<weight[1][k][i]<<endl;				
				}
			}
		}
		

		cout<<"\n Reverse plotting the traces:"<<endl;
		
		cout<<" For layers 3 and 1"<<endl;
		for(int k=0; k<58; k++)
		{
		
			if(detection[3][k]>0)
			{
				for (int i=0; i<serno[1][k]; i++)
				{
					if(mostweight[1][k]<weight[1][k][i])
					{
						mostweight[1][k]=weight[1][k][i];
						sureangle[1][k]=anglearray[1][k][i];
						surecounter[1][k]=probabletrace[1][k][i];
						surexa[1][k]=xa[1][k][i];
						sureya[1][k]=ya[1][k][i];
					}
						
				}
				cout<<"mostweight for "<<k<<" is "<<mostweight[1][k]<<endl;
				if (mostweight[1][k]==0)
				{
					if(hitsonlower[1][k]==1)
					{	
						sureangle[1][k]=anglearray[1][k][0];
						surecounter[1][k]=probabletrace[1][k][0];
						surexa[1][k]=xa[1][k][0];
						sureya[1][k]=ya[1][k][0];
						goto display2;
					}
					else
					cout<<"There isnt a deciding value"<<endl;
					
				}
				else
				{
					display2:
					cout<<"Trace between: counter# "<<k<<" and counter # "<<surecounter[1][k]<<" has angle= "<<sureangle[1][k]<<endl;
					TLine *line3 = new TLine();
					line3->SetLineWidth(2);
					line3->SetLineColor(8);
					line3->Draw();
					line3->DrawLine(surexa[1][k],sureya[1][k],xb[3][k],yb[3][k]);
				}
				
			}
			

			
					
		}		


		

		//giving the condition to draw the plot only when there are traces available
		
		for(int j=0; j<aa.np; j++)
		{
			if(trace1[j]==2 || trace2[j]==2)		
					c1[shower]->Write();
					//c1[shower]->Print("plots.pdf","Title:One bin filled");
		}


	}

		
}


