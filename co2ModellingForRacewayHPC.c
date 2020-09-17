#include "udf.h"
#include <math.h>
#include "cxndsearch.h"
#include "surf.h"

double PH = 8.0;

/* int interiorIDs[] = {3};  */
int interiorIDs[] = { 3, 10 };
int interiorIDLength;
int count;
int numberOfCells;
double phaseLevel = 0.08;
/* double phaseLevel = 0.06;  */

double lightIntensity[] = { 0,250,375,500,625,750,1000,1250,1500,2000};
/* O2 evuolution rate umol/mg Chla /h*/
double growthRate[] = {0,210,290,350,400,420,450,460,460,460};
double PHData[] = { 0, 1, 2, 3, 4, 4.500, 5, 5.500, 6, 6.370, 6.500, 7, 7.500, 7.600, 8, 8.380, 8.500, 9, 9.050, 9.500, 10, 10.30, 10.50, 11, 11.50, 12, 13, 14 };
double HCO3mData[] = { 0, 0, 0, 0, 0, 0.02000, 0.05000, 0.1350, 0.3000, 0.5000, 0.6000, 0.8100, 0.9300, 0.9500, 0.9700, 0.9760, 0.9700, 0.9630, 0.9620, 0.8900, 0.6700, 0.5000, 0.4000, 0.1800, 0.05000, 0, 0, 0 };
double CO2Data[] = { 1, 1, 1, 1, 1, 0.9800, 0.9500, 0.8650, 0.7000, 0.5000, 0.4000, 0.1900, 0.07000, 0.05000, 0.02000, 0.01200, 0.01000, 0.002000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
double CO32mData[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.01000, 0.01200, 0.02000, 0.03500, 0.03800, 0.1100, 0.3300, 0.5000, 0.6000, 0.8200, 0.9500, 1, 1, 1 };
double PHforEfficiency[]={0,6.500,6.550,6.600,6.650,6.700,6.750,6.800,6.850,6.900,6.950,7,7.050,7.100,7.150,7.200,7.250,7.300,7.350,7.400,7.450,7.500,7.550,7.600,7.650,7.700,7.750,7.800,7.850,7.900,7.950,8,8.050,8.100,8.150,8.200,8.250,8.300,8.350,8.400,8.450,8.500,100};
double Pefficiency[]={0,0,1.25737682214811e-05,3.99593527672633e-05,0.000119296591353013,0.000334575564412215,0.000881489205918614,0.00218170673761439,0.00507262014324939,0.0110796210298451,0.0227339062539777,0.0438207512339213,0.0793491295891684,0.134977416282970,0.215693297066280,0.323793989164730,0.456622713472555,0.604926811297858,0.752843580387010,0.880163316910750,0.966670292007123,0.997355701003582,0.966670292007123,0.880163316910750,0.752843580387010,0.604926811297858,0.456622713472555,0.323793989164730,0.215693297066280,0.134977416282970,0.0793491295891684,0.0438207512339213,0.0227339062539774,0.0110796210298451,0.00507262014324939,0.00218170673761443,0.000881489205918614,0.000334575564412209,0.000119296591353013,3.99593527672633e-05,1.25737682214813e-05,0,0};


int NSampleLight = sizeof(lightIntensity) / sizeof(lightIntensity[0]);
int NSamplePH = sizeof(PHData) / sizeof(PHData[0]);
int NSamplePefficiency = sizeof(PHforEfficiency) / sizeof(PHforEfficiency[0]);

double cellDensity = 60; /* mg/L */
double maxLightIntensity = 500;
static ND_Search *domain_table = NULL;

/* 
x=linspace(0,0.06,100);
 y=1./exp(250*x)*500;
 plot(x,y) 
 */

double linearIntepolationPefficiency(double PHinPut)
{
	int i;
	for (i = 0; i < NSamplePefficiency - 1; i++)
	{
		if (PHinPut < PHforEfficiency[i + 1])
			break;
	}

	if (i == NSampleLight - 1)
		i = NSampleLight - 2;

	double r = Pefficiency[i] + (Pefficiency[i] - Pefficiency[i + 1]) / (PHforEfficiency[i] - PHforEfficiency[i + 1]) *(PHinPut - PHforEfficiency[i]);
    /* printf("r=%f and  %f,  %f\n",r, Pefficiency[i],Pefficiency[i+1]); */

	return r;
}

  
double linearIntepolationLight(double intensity)
{
	int i;
	for (i = 0; i < NSampleLight - 1; i++)
	{
		if (intensity < lightIntensity[i + 1])
			break;
	}

	if (i == NSampleLight - 1)
		i = NSampleLight - 2;

	double r = growthRate[i] + (growthRate[i] - growthRate[i + 1]) / (lightIntensity[i] - lightIntensity[i + 1]) *(intensity - lightIntensity[i]);
    /* printf("r=%f and  %f,  %f\n",r, growthRate[i],growthRate[i+1]); */
	r=r/3.8333e-07* 5.0505e-08;
	return r;
	/* Normalized to Douchi, Damien’s prediction because Peer's prediction is based on O2 evolution and is not clear how to convert to C consumption.*/
} 


/*return an array[HCO3m,CO2,CO32m] */

void linearIntepolationPH(double PH, double *result)
{
	int i;
	for (i = 0; i < NSamplePH - 1; i++)
	{
		if (PH < PHData[i + 1])
			break;
	}

	if (i == NSamplePH - 1)
		i = NSamplePH - 2;

	result[0] = HCO3mData[i] + (HCO3mData[i] - HCO3mData[i + 1]) / (PHData[i] - PHData[i + 1]) *(PH - PHData[i]);
	result[1] = CO2Data[i] + (CO2Data[i] - CO2Data[i + 1]) / (PHData[i] - PHData[i + 1]) *(PH - PHData[i]);
	result[2] = 1 - result[0] - result[1];
}

DEFINE_EXECUTE_ON_LOADING(Initialize, libname)
{
	/*1. calculate how many cells are in the filter each node. 
	2. Allocate memory for filter_cell_status and current_add. 
	3. initialize filter_cell_status. 
	4. set seed for random number generator */
	interiorIDLength = sizeof(interiorIDs) / sizeof(interiorIDs[0]);
	#if RP_HOST
		printf("NSampleLight=%d\n", NSampleLight);
		printf("NSamplePH=%d\n", NSamplePH);
		printf("interiorIDLength=%d\n", interiorIDLength);
		printf("/*This is Version 1.6 */\n");
	#endif
	fflush(stdout);
}

DEFINE_EXECUTE_AT_END(reEquilibriumandSpeciesFraction)
{
	/* At the end of each time step, we need to perform:
	1.get the mass fraction of all the species.
	2. predict the reaction rate and manually predict the H consumption.
	3. get the PH at the beginning of this time step by accessing UDM.
	4. get the OH at the beginning of this time step and with the [H] at the end of this time step.
	5. build a parabola equation and pick the correct root and update the PH.
	6. reset the species fraction by the new PH.

	*/
	Domain *domain = Get_Domain(3); /*1 is for single phase */
	Thread * thread;
	cell_t c;
	count = 0;
	double hCO3, h20, CO2, CO3m2, H,OH,totalCarbon,lastPH,lastOH,lastH,bb,cc,delta,x1,x2;
	double reactionRate,reactionRateTemp,lightDepth,I;
	real x[ND_ND];
	int i,print;
#if !RP_NODE
	printf("UDF is running \n");
#endif

	for (i = 0; i < interiorIDLength; i++)
	{
		thread = Lookup_Thread(domain, interiorIDs[i]);
		begin_c_loop(c, thread)
		{
			double cVOF = C_VOF(c, thread);

			if (count<200 && x[2]>0.04&&cVOF > 0.1)
				print=1;
			else
				print=0;
			/*Those are mole concentration=massFraction/molecularMass*1000 */
			H = C_YI(c, thread, 0) / 1 * 1000;
			OH = C_YI(c, thread, 1) / 17 * 1000;
			hCO3 = C_YI(c, thread, 2) / 63 * 1000;
			CO3m2 = C_YI(c, thread, 3) / 62 * 1000;
			CO2 = C_YI(c, thread, 4) / 44 * 1000;
			h20 = C_YI(c, thread, 5) / 18 * 1000;
			totalCarbon = hCO3 + CO2 + CO3m2;
		/*get the PH at this time, It is assumed that PH does not change much during Dt. */		
			if (H<1e-17)
				PH=7;  /*if it is at air phase and H is zero, assume the PH is 7*/
			else
				PH = -log10(H); /*convert mole concentration of proton to PH */
			/* predict the proton consumption rate.*/
			double Pefficiency=(-1.61+0.47*PH-0.03*PH*PH)/0.2308333333; /*from https://link.springer.com/article/10.1007/s10811-013-0177-2 */
			if (Pefficiency<0.)
				Pefficiency=0.0;
			if (cVOF > 0.1)
			{

				C_CENTROID(x, c, thread);
				lightDepth = phaseLevel - x[2];
				if (lightDepth < 0)
					lightDepth = 0.;

				I = 1. / exp(250 *(lightDepth)) *1500;
				reactionRateTemp = linearIntepolationLight(I);
				fflush(stdout);
				reactionRate = reactionRateTemp*Pefficiency; /*mol/L/s to SI units kmol/m3/s*/
				C_UDMI(c,thread,0)=C_UDMI(c,thread,0)+reactionRate*CURRENT_TIMESTEP;
				if (hCO3>reactionRate)
					hCO3=hCO3-reactionRate*CURRENT_TIMESTEP;
				C_YI(c, thread, 2)=hCO3/1000*63; /* convert concentration to mass fraction */
 					
			}
			else
			{
				reactionRate=0.;
			}

		    if (print)
			{
				printf("count: %d my x is %f,%f,%f, cVOF=%E, lightDepth is %f, I is %E, reactionRate is %E, PH is %f, Pefficiency is %f, before consumption H is %E, OH is %E",
					count,x[0],x[1],x[2],cVOF,lightDepth,I,reactionRate,PH,Pefficiency,H,OH);
			}

			H=H-CURRENT_TIMESTEP*reactionRate; /*New [H] after proton consumption*/ 
			if (H<0)
				H=0.;

			 /*check */
		    if (print)
			{
				printf("after consumption: [H] is %E, [OH] is %E-----",H,OH);
			}
			fflush(stdout); 

		    /* solve for parabola equation*/
 			bb=H+OH;
 			cc=H*OH-1e-14;
 			delta=bb*bb-4*1*cc;
 			x1=(-bb+sqrt(delta))/2;
 			x2=(-bb-sqrt(delta))/2;
 			double testH=H+x1;
			double testOH=OH+x1;
			if (testH<0 || testOH<0)
			{
				H=H+x2;
				OH=OH+x2;
			}
			else
			{
				H=H+x1;
				OH=OH+x1;
			}
			/*check */
		    if (print)
			{
				printf("after reequrilibrum:x1=%E,x2=%E,H is %E, myid=%d  \n",x1,x2,H,myid);
				count++;
			}
			fflush(stdout);



			/* update PH */
			if (H<1e-17)
				PH=7;  /*if it is at air phase and H is zero, assume the PH is 7*/
			else
				PH = -log10(H); /*convert mole concentration of proton to PH */

			double carbonFraction[3];
			linearIntepolationPH(PH, carbonFraction);
			/*return an array[HCO3m,CO2,CO32m] */
			hCO3 = totalCarbon *carbonFraction[0];
			CO2 = totalCarbon *carbonFraction[1];
			CO3m2 = totalCarbon - hCO3 - CO2;

			h20 = (1000 - (hCO3 * 63 + CO3m2 * 62 + CO2 * 44 + H *1)) / 18.0;
			if (cVOF > 1e-8)
			{
				C_YI(c, thread, 0)=H/1000*1;
				C_YI(c, thread, 1)=OH/1000*17;
				C_YI(c, thread, 2)=hCO3/1000*63;
				C_YI(c, thread, 3)=CO3m2/1000*62;
				C_YI(c, thread, 4)=CO2/1000*44;
				C_YI(c, thread, 5) = 1 - C_YI(c, thread, 0) - C_YI(c, thread, 1) - C_YI(c, thread, 2) - C_YI(c, thread, 3) - C_YI(c, thread, 4);
			}
		}
		end_c_loop(c, thread);
	}

/*reset VOF for rotate zone*/

	thread = Lookup_Thread(domain, 10);
		begin_c_loop(c, thread)
		{		
			double cVOF = C_VOF(c, thread);
			if (cVOF<0.001)
				C_VOF(c, thread)=0.;
		}
		end_c_loop(c, thread);

	fflush(stdout);
}

DEFINE_ON_DEMAND(manuallyPatch)
{
	cell_t c;
	Domain *domain = Get_Domain(3); /*1 is for single phase */
	Thread * thread;
	printf ("manuallyPatch is running and interiorIDLength is %d\n",interiorIDLength);
	int i;
	for (i = 0; i < interiorIDLength; i++)
	{
		thread = Lookup_Thread(domain, interiorIDs[i]);
		begin_c_loop(c, thread)
		{
			double cVOF = C_VOF(c, thread);
			if (cVOF > 1e-5)
			{
				C_YI(c, thread, 0) = 1e-11;
				double H = C_YI(c, thread, 0) / 1 * 1000;
				double OH=1e-14/H;
				C_YI(c, thread, 1) = OH/1000*17;
				C_YI(c, thread, 2) = 6.3e-4;
				C_YI(c, thread, 3) = 6.2e-4;
				C_YI(c, thread, 4) = 4.4e-4;
				C_YI(c, thread, 5) = 1 - C_YI(c, thread, 0) - C_YI(c, thread, 1) - C_YI(c, thread, 2) - C_YI(c, thread, 3)- C_YI(c, thread, 4);
			}
		}
		end_c_loop(c, thread);
	}
	fflush(stdout);
}

DEFINE_ON_DEMAND(SetInitialUDM)
{
	Domain *domain = Get_Domain(3); /*1 is for single phase */
	Thread * thread;
	cell_t c;
 	double H;
	int i; 
	i=0;
 

	for (i = 0; i < interiorIDLength; i++)
	{
		thread = Lookup_Thread(domain, interiorIDs[i]);
			begin_c_loop(c, thread)
			{		
				H = C_YI(c, thread, 0) / 1 * 1000;
				if (H<1e-17)
					PH=7; /* if it is at air phase and H is zero, assume the PH is 7*/
				else
					PH = -log10(H); /*convert mole concentration of proton to PH */
				C_UDMI(c, thread, 0) = PH;
			}
			end_c_loop(c, thread);
	}
}     

/* DEFINE_ON_DEMAND(resetDensity)
{
	cell_t c;

	double hCO3, h20, CO2, CO3m2, H, totalCarbon;
	Domain *domain = Get_Domain(3); 
	Thread * thread;
	printf ("manuallyPatch is running\n");
		thread = Lookup_Thread(domain, 16);
		begin_c_loop(c, thread)
		{		
			double cVOF = C_VOF(c, thread);
			if (cVOF<0.001)
				C_VOF(c, thread)=0.;


		}

		end_c_loop(c, thread);
	fflush(stdout);
}     */



DEFINE_ON_DEMAND(testlinearIntepolationPH11111)
{
	double result[10];
	linearIntepolationPH(8.0, result);
	 printf("%f;%f;%f\n",result[0],result[1],result[2]);
	fflush(stdout);
}   


DEFINE_ON_DEMAND(resetSpeciesFractionandReEquilibriumByPH)
{
	/* At the end of each time step, we need to perform:
	1.get the mass fraction of all the species.
	2. predict the reaction rate and manually predict the H consumption.
	3. get the PH at the beginning of this time step by accessing UDM.
	4. get the OH at the beginning of this time step and with the [H] at the end of this time step.
	5. build a parabola equation and pick the correct root and update the PH.
	6. reset the species fraction by the new PH.

	*/
}

 

DEFINE_LINEARIZED_MASS_TRANSFER(co2Loss, cell, thread, from_index, from_species_index, to_index, to_species_index, d_mdot_d_vof_from, d_mdot_d_vof_to)
{
	Thread *liq = THREAD_SUB_THREAD(thread, from_index);
	double cVOF, massTransferRate;
	double CO2;

	massTransferRate = 0.0;
	cVOF = C_VOF(cell, liq);

	if (cVOF >= 0.001 && cVOF <= 0.9)
	{
		CO2 = C_YI(cell, liq, 3) / 44 * 1000;
		massTransferRate = CO2 * 10;
		fflush(stdout);
	}
	return (massTransferRate);
}


DEFINE_HET_RXN_RATE(consumption, c, t, hr, mw, yi, rr, rr_t)
{
	/* excuted on each cell at each iteration
	2. Allocate memory for filter_cell_status and current_add. 
	3. initialize filter_cell_status. 
	4. set seed for random number generator */

	Thread **pt = THREAD_SUB_THREADS(t);
	Thread *tp = pt[0];	 /* primary */
	Thread *ts = pt[1];	 /* secondary */
	double cVOF = C_VOF(c, ts);
	double H=C_YI(c, ts, 0) / 1 * 1000;
	if (H<1e-17)
		PH=7.123;  /*if it is at air phase and H is zero, assume the PH is 7*/
	else
		PH = -log10(H); /*convert mole concentration of proton to PH */
 	/* double Pefficiency=linearIntepolationPefficiency(PH); */
	double Pefficiency=(-1.61+0.47*PH-0.03*PH*PH)/0.2308333333; /*from https://link.springer.com/article/10.1007/s10811-013-0177-2 */
	if (Pefficiency<0.)
		Pefficiency=0.0;

	/* printf("I am cell %d, my ph is %0.6e and my Pefficiency is %0.6e\n",c,PH,Pefficiency); */
	if (cVOF > 0.1)
	{
		real x[ND_ND];
		C_CENTROID(x, c, t);
		double lightDepth = phaseLevel - x[2];
		if (lightDepth < 0)
			lightDepth = 0.;

		double I = 1. / exp(250 *(lightDepth)) *1500;
		double reactionRate = linearIntepolationLight(I);
		fflush(stdout);
		*rr = reactionRate/1000000*Pefficiency; /*mol/L/s to SI units kmol/m3/s*/
		*rr=0; /* disable reaction, do it manually in DEFINE_EXECUTE_AT_END(reEquilibriumandSpeciesFraction)
		/*
		if (count<10)
		{
			printf("my x is %f, lightDepth is %f, I is %f, rris %f, PH is %f, Pefficiency is %f  \n",x[1],lightDepth,I,*rr,PH,Pefficiency);
			count++;
		}
		fflush(stdout);
		*/
	}
	else
	{
		*rr = 0.;
	}
}

DEFINE_LINEARIZED_MASS_TRANSFER(evaporation, cell, thread, from_index, from_species_index, to_index, to_species_index, d_mdot_d_vof_from, d_mdot_d_vof_to)
{
	Thread *liq = THREAD_SUB_THREAD(thread, from_index);
	Thread *gas = THREAD_SUB_THREAD(thread, to_index);
	int cID;
	double cVOF, massTransferRate;
	
	massTransferRate = 0.0;
	cID = cell;
	cVOF = C_VOF(cell, liq);
	if (cVOF < 0.9)
	{
		massTransferRate = 0;
	}

	fflush(stdout);
	return (massTransferRate);
}
 