#include "udf.h"
#include <math.h>
#include "cxndsearch.h"
#include "surf.h"

double PH = 8.0;

/* int interiorIDs[] = {3};  */
int interiorIDs[] = { 3, 12, 16 };
int interiorIDLength;
int count;
int numberOfCells;
double phaseLevel = 0.20;
/* double phaseLevel = 0.06;  */

double lightIntensity[] = { 0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, 400, 500, 800, 1000, 2000 };
double growthRate[] = { 0, 2.50000000000000, 2.70000000000000, 2, 1.50000000000000, 1.30000000000000, 1.20000000000000, 1.10000000000000, 1, 0.900000000000000, 0.800000000000000, 0.750000000000000, 0.730000000000000, 0.710000000000000, 0.680000000000000, 0.600000000000000, 0.550000000000000, 0.500000000000000, 0.400000000000000, 0.300000000000000, 0.300000000000000 };
double PHData[] = { 0, 1, 2, 3, 4, 4.50000000000000, 5, 5.50000000000000, 6, 6.37000000000000, 6.50000000000000, 7, 7.50000000000000, 7.60000000000000, 8, 8.38000000000000, 8.50000000000000, 9, 9.05000000000000, 9.50000000000000, 10, 10.3000000000000, 10.5000000000000, 11, 11.5000000000000, 12, 13, 14 };
double HCO3mData[] = { 0, 0, 0, 0, 0, 0.0200000000000000, 0.0500000000000000, 0.135000000000000, 0.300000000000000, 0.500000000000000, 0.600000000000000, 0.810000000000000, 0.930000000000000, 0.950000000000000, 0.970000000000000, 0.976000000000000, 0.970000000000000, 0.963000000000000, 0.962000000000000, 0.890000000000000, 0.670000000000000, 0.500000000000000, 0.400000000000000, 0.180000000000000, 0.0500000000000000, 0, 0, 0 };
double CO2Data[] = { 1, 1, 1, 1, 1, 0.980000000000000, 0.950000000000000, 0.865000000000000, 0.700000000000000, 0.500000000000000, 0.400000000000000, 0.190000000000000, 0.0700000000000000, 0.0500000000000000, 0.0200000000000000, 0.0120000000000000, 0.0100000000000000, 0.00200000000000000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
double CO32mData[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0100000000000000, 0.0120000000000000, 0.0200000000000000, 0.0350000000000000, 0.0380000000000000, 0.110000000000000, 0.330000000000000, 0.500000000000000, 0.600000000000000, 0.820000000000000, 0.950000000000000, 1, 1, 1 };

int NSampleLight = sizeof(lightIntensity) / sizeof(lightIntensity[0]);
int NSamplePH = sizeof(PHData) / sizeof(PHData[0]);

double maxLightIntensity = 500;
static ND_Search *domain_table = NULL;

/* 
x=linspace(0,0.06,100);
 y=1./exp(250*x)*500;
 plot(x,y) 
 */

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

	return r;
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
		printf("/*This is Version 1.4 */\n");
	#endif
	fflush(stdout);
}

DEFINE_EXECUTE_AT_END(resetSpeciesFraction)
{
	Domain *domain = Get_Domain(3); /*1 is for single phase */
	Thread * thread;
	cell_t c;
	count = 0;
	double hCO3, h20, CO2, CO3m2, H, totalCarbon;
	int i;

	for (i = 0; i < interiorIDLength; i++)
	{
		thread = Lookup_Thread(domain, interiorIDs[i]);
		begin_c_loop(c, thread)
		{
			double cVOF = C_VOF(c, thread);
			/*Those are mole concentration=massFraction/molecularMass*1000 */
			H = C_YI(c, thread, 0) / 1 * 1000;
			hCO3 = C_YI(c, thread, 1) / 63 * 1000;
			CO3m2 = C_YI(c, thread, 2) / 62 * 1000;
			CO2 = C_YI(c, thread, 3) / 44 * 1000;
			h20 = C_YI(c, thread, 4) / 18 * 1000;
			totalCarbon = hCO3 + CO2 + CO3m2;
			/*equilibrium mole fraction*/
			if (H<1e-17)
				PH=7; /* if it is at air phase and H is zero, assume the PH is 7*/
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
				C_YI(c, thread, 1)=hCO3/1000*63;
				C_YI(c, thread, 2)=CO3m2/1000*62;
				C_YI(c, thread, 3)=CO2/1000*44;
				C_YI(c, thread, 4) = 1 - C_YI(c, thread, 0) - C_YI(c, thread, 1) - C_YI(c, thread, 2) - C_YI(c, thread, 3);
			}
		}

		end_c_loop(c, thread);
	}

	fflush(stdout);
}

DEFINE_ON_DEMAND(manuallyPatch)
{
	cell_t c;
	Domain *domain = Get_Domain(3); /*1 is for single phase */
	Thread * thread;
	printf ("manuallyPatch is running\n");
	int i;
	for (i = 0; i < interiorIDLength; i++)
	{
		thread = Lookup_Thread(domain, interiorIDs[i]);
		begin_c_loop(c, thread)
		{
			double cVOF = C_VOF(c, thread);
			if (cVOF > 1e-8)
			{
				C_YI(c, thread, 0) = 1e-13;
				C_YI(c, thread, 1) = 0.001;
				C_YI(c, thread, 2) = 0.002;
				C_YI(c, thread, 3) = 0.003;
				C_YI(c, thread, 4) = 1 - C_YI(c, thread, 0) - C_YI(c, thread, 1) - C_YI(c, thread, 2) - C_YI(c, thread, 3);
			}
		}
		end_c_loop(c, thread);
	}
	fflush(stdout);
}


/*   DEFINE_ON_DEMAND(resetDensity)
{
	cell_t c;
	double hCO3, h20, CO2, CO3m2, H, totalCarbon;
	Domain *domain = Get_Domain(1);  
	Thread * thread;
	printf ("manuallyPatch is running\n");
	for (int i = 0; i < interiorIDLength; i++)
	{
		thread = Lookup_Thread(domain, interiorIDs[i]);
		begin_c_loop(c, thread)
		{
				C_R(c, thread) = 1998.22233;
		}

		end_c_loop(c, thread);
	}
	fflush(stdout);
}  

			*/



DEFINE_ON_DEMAND(resetSpeciesFractionByPH)
{
	Domain *domain = Get_Domain(3); /*1 is for single phase */
	Thread * thread;
	cell_t c;
	double hCO3, h20, CO2, CO3m2, H, totalCarbon;
	int i;
	for (i = 0; i < interiorIDLength; i++)
	{
		thread = Lookup_Thread(domain, interiorIDs[i]);
		begin_c_loop(c, thread)
		{
			double cVOF = C_VOF(c, thread);

			/*Those are mole concentration=massFraction/molecularMass*1000 */
			H = C_YI(c, thread, 0) / 1 * 1000;
			hCO3 = C_YI(c, thread, 1) / 63 * 1000;
			CO3m2 = C_YI(c, thread, 2) / 62 * 1000;
			CO2 = C_YI(c, thread, 3) / 44 * 1000;
			h20 = C_YI(c, thread, 4) / 18 * 1000;
			totalCarbon = hCO3 + CO2 + CO3m2;

			/*equilibrium mole fraction*/
			if (H<1e-17)
				PH=7; /* if it is at air phase and H is zero, assume the PH is 7*/
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
				C_YI(c, thread, 1)=hCO3/1000*63;
				C_YI(c, thread, 2)=CO3m2/1000*62;
				C_YI(c, thread, 3)=CO2/1000*44;
				C_YI(c, thread, 4) = 1 - C_YI(c, thread, 0) - C_YI(c, thread, 1) - C_YI(c, thread, 2) - C_YI(c, thread, 3);
			}
		}

		end_c_loop(c, thread);
	}
}

DEFINE_ON_DEMAND(resetCount)
{
	count = 0;
	printf("count is reset as %d \n", count);
	fflush(stdout);

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

	if (cVOF > 0.1)
	{
		real x[ND_ND];
		C_CENTROID(x, c, t);
		double lightDepth = phaseLevel - x[2];
		if (lightDepth < 0)
			lightDepth = 0.;

		double I = 1. / exp(250 *(lightDepth)) *500;
		double reactionRate = linearIntepolationLight(I);
		fflush(stdout);
		*rr = reactionRate;
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

 