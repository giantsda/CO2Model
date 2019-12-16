#include "udf.h"
#include <math.h>
#include "cxndsearch.h"
#include "surf.h"


#define K1 2.0e-2
#define K2 5.

double PH=8.0;
double* totalCarbon;
int interiorID=3;
int numberOfCells;
double phaseLevel=0.06;
int haha=0;
int count=0;
int counttt=0;

double lightIntensity[]= {0,25,50,75,100,125,150,175,200,225,250,275,300,325,350,375,400,500,800,1000,2000};
double growthRate[]= {0,2.50000000000000,2.70000000000000,2,1.50000000000000,1.30000000000000,1.20000000000000,1.10000000000000,1,0.900000000000000,0.800000000000000,0.750000000000000,0.730000000000000,0.710000000000000,0.680000000000000,0.600000000000000,0.550000000000000,0.500000000000000,0.400000000000000,0.300000000000000,0.300000000000000};
int NSample=sizeof(lightIntensity)/sizeof(lightIntensity[0]);

double maxLightIntensity=500;
static ND_Search *domain_table=NULL;


/*  
x=linspace(0,0.06,100);
 y=1./exp(250*x)*500;
 plot(x,y) 
 */

double linearIntepolation(double intensity)
	{
 
		int i;
		for(i=0;i<NSample-1;i++) 
			{
				if (intensity<lightIntensity[i+1])
					break;
			}
		if (i==NSample-1)	
			i=NSample-2;

		// printf("i=%d\n",i); 
		double r=growthRate[i]+(growthRate[i]-growthRate[i+1])/(lightIntensity[i]-lightIntensity[i+1])*(intensity-lightIntensity[i]);
		// printf("r=%f and  %f ,  %f\n",r, growthRate[i],growthRate[i+1]); 

		return r;
	}



double getHCO3FractionByPH(double PH)
{
		return 1./3;

}

double getCO32mFractionByPH(double PH)
{

		return 1./3;
}

DEFINE_EXECUTE_ON_LOADING(Initialize, libname) 
{
	/* 1. calculate how many cells are in the filter each node. 
	2. Allocate memory for filter_cell_status and current_add. 
	3. initialize filter_cell_status. 
	4. set seed for random number generator */
  numberOfCells = 0;
  Domain * domain = Get_Domain (1);  /* 1 is for single phase */
  Thread * thread = Lookup_Thread (domain, interiorID); 
  cell_t c;
  // printf ("Initialize si running\n");
  // fflush (stdout);

  begin_c_loop (c, thread)
    {
      numberOfCells++;
    }
  end_c_loop (c, thread);
  printf ("\n numberOfCells=%d on node:%d <<<<<<<<<<<<<<<<<<<< \n",
	  numberOfCells, myid);
  totalCarbon = (double*) malloc (sizeof(double) * numberOfCells);
  // current_add = (double*) malloc (sizeof(double) * numberOfCells);
  int i;
  for (i = 0; i < numberOfCells; i++) /* initialize totalCarbon */
    {
      totalCarbon[i] = 0.;
    }
  // srand (time (0) + myid * 123);  /* set seed for random number */
	printf("NSample=%d\n",NSample);
	fflush (stdout);
}


// DEFINE_PROFILE(species_profile, t, i)
// {
// cell_t c;

// begin_c_loop(c, t)
// {
// 		{
// 			C_PROFILE(c,t,i) = 0.2;
// 		}
// }
// end_c_loop(c, t)
// }

 

DEFINE_EXECUTE_AT_END(resetSpeciesFraction)
{
	Domain * domain = Get_Domain (3);  /* 1 is for single phase */
	Thread * thread = Lookup_Thread (domain, interiorID); 
	cell_t c;
// #if RP_NODE
	double hCO3,h20,CO2,CO3m2,totalCarbon;

  // printf ("resetSpeciesFraction si running\n");
  // fflush (stdout);


	begin_c_loop(c, thread) 
		{	
			int i=0;
			/* Those are mole concentration=massFraction/molecularMass*1000 */
				 hCO3=C_YI(c, thread, 0)/63*1000;
				 h20=C_YI(c, thread, 1)/18*1000;
				 CO2=C_YI(c, thread, 2)/44*1000;
				 CO3m2=C_YI(c, thread, 3)/62*1000;
				 totalCarbon=hCO3+CO2+CO3m2;
			/* equilibrium mole fraction*/	 
				 hCO3=totalCarbon*getHCO3FractionByPH(PH);
				 CO3m2=totalCarbon*getCO32mFractionByPH(PH);
				 CO2=totalCarbon-hCO3-CO3m2;
				 h20=(1000-(hCO3*63+CO3m2*62+CO2*44))/18.0;
		    /*convert to mass fraction=moleFraction/1000*molecularMass*/
				C_YI(c, thread, 0)=hCO3/1000*63;
				C_YI(c, thread, 1)=h20/1000*18;
				C_YI(c, thread, 2)=CO2/1000*44;
				C_YI(c, thread, 3)=CO3m2/1000*62;
 
			// if (c<=100)
			// {
			// 	// printf("hCO3 is %f for cell %d in node %d \n",hCO3,c,myid);
			// 	// printf("h20 is %f for cell %d in node %d \n",h20,c,myid);
			// 	// printf("CO2 is %f for cell %d in node %d \n",CO2,c,myid);
			// 	// printf("CO3m2 is %f for cell %d in node %d \n",CO3m2,c,myid);
			// 	printf("hCO3 is %f for cell %d in node %d \n",C_YI(c, thread, 0),c,myid);
			// 	printf("h20 is %f for cell %d in node %d \n",C_YI(c, thread, 1),c,myid);
			// 	printf("CO2 is %f for cell %d in node %d \n",C_YI(c, thread, 2),c,myid);
			// 	printf("CO3m2 is %f for cell %d in node %d \n",C_YI(c, thread, 3),c,myid);
			// }
		} 

	end_c_loop (c, thread);
 
// #endif
	// printf("thread is %p \n",thread);
	fflush(stdout);
}

DEFINE_ON_DEMAND(resetSpeciesFractionByPH)
{
	Domain * domain = Get_Domain (3);  /* 1 is for single phase */
	Thread * thread = Lookup_Thread (domain, interiorID); 
	cell_t c;
// #if RP_NODE
	double hCO3,h20,CO2,CO3m2,totalCarbon;

  // printf ("resetSpeciesFractionByPH si running\n");
  // fflush (stdout);
	
	// Thread **pt = THREAD_SUB_THREADS(thread);
	// Thread *tp = pt[0]; // primary	
	// Thread *ts = pt[1]; // secondary
 
	begin_c_loop(c, thread) 
		{	
				double cVOF=C_VOF(c,thread);

			/* Those are mole concentration=massFraction/molecularMass*1000 */
				 hCO3=C_YI(c, thread, 0)/63*1000;
				 h20=C_YI(c, thread, 1)/18*1000;
				 CO2=C_YI(c, thread, 2)/44*1000;
				 CO3m2=C_YI(c, thread, 3)/62*1000;
				 totalCarbon=hCO3+CO2+CO3m2;
			/* equilibrium mole fraction*/	 
				 hCO3=totalCarbon*getHCO3FractionByPH(PH);
				 CO3m2=totalCarbon*getCO32mFractionByPH(PH);
				 CO2=totalCarbon-hCO3-CO3m2;
				 h20=(1000-(hCO3*63+CO3m2*62+CO2*44))/18.0;
		    /*convert to mass fraction=moleFraction/1000*molecularMass*/
				C_YI(c, thread, 0)=hCO3/1000*63;
				C_YI(c, thread, 1)=h20/1000*18;
				C_YI(c, thread, 2)=CO2/1000*44;
				C_YI(c, thread, 3)=CO3m2/1000*62;
 
			if (c<=100)
			{
				// printf("hCO3 is %f for cell %d in node %d \n",hCO3,c,myid);
				// printf("h20 is %f for cell %d in node %d \n",h20,c,myid);
				// printf("CO2 is %f for cell %d in node %d \n",CO2,c,myid);
				// printf("CO3m2 is %f for cell %d in node %d \n",CO3m2,c,myid);
				printf("hCO3 is %f for cell %d , cVOF=%f, in node %d \n",C_YI(c, thread, 0),c,cVOF,myid);
				printf("h20 is %f for cell %d in node %d \n",C_YI(c, thread, 1),c,myid);
				printf("CO2 is %f for cell %d in node %d \n",C_YI(c, thread, 2),c,myid);
				printf("CO3m2 is %f for cell %d in node %d \n",C_YI(c, thread, 3),c,myid);
			}
		}
	end_c_loop (c, thread);
}

DEFINE_ON_DEMAND(testFindCell)
{

  printf ("testFindCell si running\n");
  fflush (stdout);

	CX_Cell_Id *cx_cell;
	cell_t c;
	Thread *t;
	real P[3];
 
	real P_Cell[3];

	P[0]=0.1;
	P[1]=0.1;
	P[2]=0.1;


	P_Cell[0]=0.1;
	P_Cell[1]=0.1;
	P_Cell[2]=0.1;


	printf("Exeed >>>>>>>> in node %d \n",myid);
	domain_table = CX_Start_ND_Point_Search( domain_table,TRUE,-1);
	float dist=0.1;
	cx_cell =CX_Find_Closest_Cell_To_Point(domain_table, P, &dist, 111111);
	printf("cx_cell is %p in node %d \n",cx_cell,myid);


	if (cx_cell) 
		{
			c = RP_CELL(cx_cell);
			t = RP_THREAD(cx_cell);
			C_CENTROID(P_Cell,c,t);
			printf("Found cell at [%g,%g,%g] with centroid [%g,%g,%g]. dist=%f, in node %d \n",P[0],P[1],P[2],P_Cell[0],P_Cell[1],P_Cell[2],dist,myid);
		} 
	else 
		{
			printf("Could not find cell at [%g,%g,%g]! in node %d \n",P[0],P[1],P[2],myid);
		}

	domain_table = CX_End_ND_Point_Search(domain_table);


	fflush (stdout);
 
}


DEFINE_LINEARIZED_MASS_TRANSFER(co2Loss, cell, thread, from_index,from_species_index, to_index, to_species_index, d_mdot_d_vof_from,d_mdot_d_vof_to)
{
   real m_lg;
   // real T_SAT = 373.15;
   Thread *liq = THREAD_SUB_THREAD(thread, from_index);
   Thread *gas = THREAD_SUB_THREAD(thread, to_index);
  
  // printf ("co2Loss is running\n");
  // fflush (stdout);

  double cVOF,massTransferRate;
  double verticalDistance=0.02;
  double xp,yp,zp;
  double CO2;
   m_lg = 0.;
  massTransferRate=0.0;
#if RP_NODE
   cVOF=C_VOF(cell,liq);
   if (cVOF>=0.1 && cVOF<=0.9)
   {
	real x[ND_ND];
	C_CENTROID(x,cell,thread); 
 	xp=x[0];
 	yp=x[1];
 	zp=x[2]-verticalDistance;

	// printf("counttt=%d \n",counttt);
	// counttt++;

	real P[3];
	real P_Cell[3];
	P[0]=xp;
	P[1]=yp;
	P[2]=zp;

	cell_t c;
	Thread *t;
	CX_Cell_Id *cx_cell;
 
	domain_table = CX_Start_ND_Point_Search( domain_table,TRUE,-1);
	float dist=0.1;
	cx_cell =CX_Find_Closest_Cell_To_Point(domain_table, P, &dist, 0.1);
	// printf("cx_cell is %p in node %d \n",cx_cell,myid);


	if (cx_cell) 
		{
			Domain * domain = Get_Domain (3);  /* 1 is for single phase */
			Thread * threadMix = Lookup_Thread (domain, interiorID); 


			c = RP_CELL(cx_cell);
			t = RP_THREAD(cx_cell);
			C_CENTROID(P_Cell,c,t);

				// if (count<=20)
				// 	{
						// printf("target location is %f,%f,%f \n",x[0],x[1],x[2]);
						// printf("c is %d and t is %p \n",c,t);
						// printf("the location by finding t is %f, %f, %f; myid=%d \n", P_Cell[0],P_Cell[1],P_Cell[2],myid);
						// printf("Found cell at [%g,%g,%g] with centroid [%g,%g,%g]. dist=%f, in node %d \n",P[0],P[1],P[2],P_Cell[0],P_Cell[1],P_Cell[2],dist,myid);
			            CO2=fabs(C_YI_G(c, threadMix,2)[2])/44*1000;
						count++;
						massTransferRate=CO2*0.01;

			            // C_YI_G(c,t)[i]


			            // C_CENTROID(P_Cell,c,threadMix);
						// printf("Co2 is %f at cell %d threadMix=%p and count is %d\n",CO2,c,threadMix,count);
						// printf("the location by threadMix is %f, %f, %f \n", P_Cell[0],P_Cell[1],P_Cell[2]);
						// printf("--------------------------------\n");

					// }
		} 
	// else 

	// 	{
	// 		printf("Could not find cell at [%g,%g,%g]! in node %d \n",P[0],P[1],P[2],myid);
	// 	}
	domain_table = CX_End_ND_Point_Search(domain_table);

 //   	massTransferRate=10.8;
   }
#endif

// if (myid==1)
// {
//       printf("this is cell %d cvof=%f from ID %d \n",cID,cVOF,myid);
// }
   
   return (massTransferRate);
}
 
DEFINE_LINEARIZED_MASS_TRANSFER(co2LossV2, cell, thread, from_index,from_species_index, to_index, to_species_index, d_mdot_d_vof_from,d_mdot_d_vof_to)
{
   real m_lg;
   // real T_SAT = 373.15;
   Thread *liq = THREAD_SUB_THREAD(thread, from_index);
	Thread *gas = THREAD_SUB_THREAD(thread, to_index);

	// printf ("co2Loss is running\n");
	// fflush (stdout);
	Domain * domain = Get_Domain (3);  /* 1 is for single phase */
	Thread * threadMix = Lookup_Thread (domain, interiorID); 
	double cVOF,massTransferRate;
	double verticalDistance=0.02;
	double xp,yp,zp;
	double CO2Gradient;
	m_lg = 0.;
    massTransferRate=0.0;
#if RP_NODE
   cVOF=C_VOF(cell,liq);
   if (cVOF>=0.1 && cVOF<=0.9)
   {
		CO2Gradient=fabs(C_YI_G(cell, threadMix,2)[2])/44*1000;
		massTransferRate=CO2*0.01;
	} 
#endif
 
   return (massTransferRate);
}


 
DEFINE_HET_RXN_RATE(consumption, c, t, hr, mw, yi, rr, rr_t)
{
	/* excuted on each cell at each iteration
	2. Allocate memory for filter_cell_status and current_add. 
	3. initialize filter_cell_status. 
	4. set seed for random number generator */

	Thread **pt = THREAD_SUB_THREADS(t);
	Thread *tp = pt[0]; // primary	
	Thread *ts = pt[1]; // secondary
  
  // printf ("consumption is running\n");
  // fflush (stdout);

	double cVOF=C_VOF(c,ts);
 
// printf("for cell %d in node %d; cVOF=%f \n",c,myid,cVOF);

	if (cVOF>0.1)
		{ 
			real x[ND_ND];
			C_CENTROID(x,c,t); 
		 	double lightDepth=phaseLevel-x[2];
		 	if (lightDepth<0)
		 		lightDepth=0.;

		 	double I=1./exp(250*(lightDepth))*500;
			// if (c<=100)
			// {
			// 	printf("z:%f for cell %d in node %d; I=%f \n",x[2],c,myid,I);
			// }
		 	double reactionRate=linearIntepolation(I);

			// if (c<=100)
			// {
			// 	printf("z:%f for cell %d in node %d; I=%f; reactionRate=%f \n",x[2],c,myid,I,reactionRate);
			// }

		 	fflush (stdout);
			*rr = reactionRate;
		}
	else
		{
			*rr=0.;
		}
}


DEFINE_LINEARIZED_MASS_TRANSFER(evaporation, cell, thread, from_index,from_species_index, to_index, to_species_index, d_mdot_d_vof_from,d_mdot_d_vof_to)
{
   real m_lg;
   // real T_SAT = 373.15;
   Thread *liq = THREAD_SUB_THREAD(thread, from_index);
   Thread *gas = THREAD_SUB_THREAD(thread, to_index);
  int cID;
  double cVOF,massTransferRate;
   // printf("this is cell  111111111111111111111111111\n");
   m_lg = 0.;

  // printf ("evaporation is running\n");
  // fflush (stdout);

	// printf("this is cell %d cvof=%f from ID %d \n",cID,cVOF,myid);

   massTransferRate=0.0;
#if RP_NODE
   cID=cell;
   cVOF=C_VOF(cell,liq);
   if (cVOF<0.9)
   {
      		// printf("this is cell %d cvof=%f from ID %d \n",cID,cVOF,myid);
   	massTransferRate=0;
   }
#endif

// if (myid==1)
// {
//       printf("this is cell %d cvof=%f from ID %d \n",cID,cVOF,myid);
// }
   
	fflush (stdout);
   return (massTransferRate);
}


