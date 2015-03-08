#include <iostream>
#include <math.h>
#include <time.h>
#include <fstream>
#include <stdio.h>
#include <string.h>
#ifdef MYRTLE
#endif

#ifdef __MPI__
/* include MPI definitions */
#include <mpi.h>
#endif

#include "macros.h"
#include "stringfunc.h"
#include "streamfunc.h"
#include "data_array.h"
#include "polygon.h"
#include "ran01.h"
#include "neuron.h"
#include "glia.h"
#include "infobjarray.h"
#include "data.h"
#include "collisions.h"
#include "defs.h"
#include "funcs.h"
#include "grid.h"
#include "move.h"
#include "pair_corr.h"

#include "generate.h"

/*
 * PROGRAM TO GENERATE A 3D NEURONAL BLOCK OF MICROCOLUMNS AND
 * THEN RANDOMIZE THEIR MOVEMENTS IN SEVERAL WAYS, IN OVERLAP AND
 * NON-OVERLAP MANNER.
 */

int main(int argc, char* argv[])
{
  int i,j,k;
  int tid, nthreads, NCPUs;
  int cf;
  int n = 4;
  int prev_time_label;


#ifdef __MPI__
  // =============================================
  // INITIALIZE MPI
    MPI_Status stat;
    
    /* add in MPI startup routines */
    /* 1st: launch the MPI processes on each node */
    MPI_Init(&argc,&argv);	

    /* 2nd: request a thread id, sometimes called a "rank" from 
            the MPI master process, which has rank or tid == 0 */
    MPI_Comm_rank(MPI_COMM_WORLD, &tid);  

    /* 3rd: this is often useful, get the number of threads
            or processes launched by MPI, this should be NCPUs-1 */
    MPI_Comm_size(MPI_COMM_WORLD, &nthreads);
	NCPUs = nthreads;

	cout << "num procs " << nthreads << "\n";
	cout << "my id     " << tid << "\n";
  // =============================================
#endif


  // CHECK FOR INPUT FILE
  if (argc<3) 
  {
	cout <<"\nThis is 'slicegen' program version 0.13\n\n";
	cout << "Usage:\n";
	cout << "\ta.out " 
	  << " <InputFile> "
	  << " <HowMany> "
	  << "\n\n";
	exit(1);
  }

  // =============================================
  // HOW MANY CONFIGURATIONS TO GENERATE. EACH ONE IS ONLY DIFFERENT
  // BY THE RANDOM GENERATOR SEED USED.
  int howmany = atoi( argv[2] );

  // =============================================
  // MAIN OBJECT HOLDING THE SIMULATION PARAMETERS
  data dq;
  // READ INPUT PARAMS
  dq.readFile(argv[1]);
  cout << "Read from data file: " << argv[1] << "\n" ;

  // =============================================
  // DEFAULT SEED FOR RANDOM GENERATOR
  long seed;                          
  //long time_sub = 1149000000;
  long time_sub = 1234234670;
  // =============================================
  // RANDOM SEED GIVEN OR GENERATED. ALLOWS FOR REPEATING RUNS EXACTLY.
  // CHECK IF PARAMETER FILE CAME WITH SEED, AND GIVE PRIORITY.
  if (dq.seed!=-1) 
  {
	seed = dq.seed;
  }
  else 
  {
  	// USE TIME AS SEED
#ifndef __MPI__
  	seed = - long( time(NULL) - time_sub );
#else
	// HAVE A DIFFERENT SEED PER PROCESSOR
  	seed = - long( time(NULL) - time_sub ) + tid;
#endif
	dq.seed = seed;
  }

  // ANOUNCE SEED IN USE
#ifndef __MPI__
  cout << "Using " << dq.seed << " as configuration seed.\n";
#else
  cout << "PROC " << tid << " ->Using " << dq.seed << " as configuration seed.\n";
#endif
  // SEED THE GENERATOR
  srand(seed);

  // =============================================
  // LOOP TO GENERATE howmany CONFIGURATIONS.
#ifdef __MPI__
  howmany /= int(nthreads);
  for (cf=0; cf<howmany; cf++)
#else
  for (cf=0; cf<howmany; cf++)
#endif
  {
#ifdef __MPI__
	cout << "*** I AM PROC " << tid << " DOING " << cf << " ***\n";
#endif
  	// READ INPUT PARAMS AGAIN FOR SUBSEQUENT
	// RUNS.
	if (cf!=0)
  		dq.readFile(argv[1]);

	// =============================================
	// GENERATE RANDOM ANGLES IF NECESSARY
	gen_rand_angles(dq);
	// INITIALIZE ARRAY OF COLLISIONS
	collisions col;
	col.init(dq.num_steps);

	// =============================================
	// GET THE TIME AND USE AS A STAMP FOR THE FILES TO BE
	// WRITTEN BELOW. SUCCESIVE INSTANCES OF THE PROGRAM THEN WILL
	// NOT OVERWRITE PREVIOUS FILES.
	time_t *the_time = new time_t[20];
	time(the_time);
	// MAKE SURE WE ARE USING A NEW LABEL
	while (*the_time == prev_time_label && cf!=0)
	{
	  cout << "waiting.............." << *the_time << " " << prev_time_label << "\n";
	  int max_iter = 5000;
	  for (int s=0; s<max_iter; s++)
		for (int t=0; t<max_iter; t++)
		{
		  float a = sqrt(sqrt(s * t * 1.0));
		}
	  time(the_time);
	}
	prev_time_label = *the_time;
#ifndef __MPI__
	// MAKE THE LABEL A LITTLE SHORTER
	int trunc_sec = *the_time - time_sub;
#else
	// USE SIMPLE RANDOM LABEL FOR FILES
	int trunc_sec = rand();
#endif
	char* file_lab = new char[10];
	itoa(trunc_sec, file_lab);
	cout << "Label identifier is: " << trunc_sec << "\n";

	// =============================================
	// WRITE INPUT PARAMS TO FILE
	char pre_par[40]     = "InputParams_";
	strcat(pre_par, file_lab);
	strcat(pre_par, "_");
	char *param_file     = NULL;
	param_file     = conc_char_2_int(pre_par,dq.num_steps,3);
	strcat(param_file,     ".txt");
	ofstream param_out(param_file);

	param_out << dq.seed  << "\tRandom seed \n";
	param_out << file_lab << "\tLabel \n";
	param_out << argv[1]  << "\tInput file \n";
	param_out << dq.inter_column << " \tinter_column\n";
	param_out << dq.put_random_col_xy << " \tput_random_col_xy\n";
	param_out << dq.random_col_mindist << " \trandom_col_mindist\n";
	param_out << dq.inter_neuron << " \tinter_neuron\n";
	param_out << dq.put_nice_z   << " \tput nice z \n";
	param_out << dq.delta_neuron << " \tdelta_neuron\n";
	param_out << dq.delta_xy_neuron << " \tdelta_xy_neuron\n";
	param_out << dq.dev_base     << " \tbase deviations \n";
	param_out << dq.skip_neurons     << " \tmissing neurons\n";
	param_out << dq.neu_radius   << " \tneu_radius  \n";
	param_out << dq.put_glia   << " \tput_glia  \n";
	param_out << dq.glia_radius   << " \tglia_radius  \n";
	param_out << dq.glia_amount   << " \tglia_amount  \n";
	param_out << dq.type_glia_position   << " \ttype_glia_position  \n";
	param_out << dq.L_tissue     << " \tL_tissue    \n";
	param_out << dq.tissue_thick << " \ttissue_thick\n";
	param_out << dq.num_steps    << " \tnum_steps   \n";
	param_out << dq.move_step    << " \tmove_step   \n";
	param_out << DEGREE(dq.angle      ) << " \tangle (blade)   \n";
	param_out << DEGREE(dq.angle_phi  ) << " \tangle_phi\n";
	param_out << dq.angle_rndphi_low  << " \tangle_rndphi_low\n";
	param_out << dq.angle_rndphi_high << " \tangle_rndphi_high\n";
	param_out << DEGREE(dq.angle_theta) << " \tangle_theta\n";
	param_out << dq.overlap      << " \toverlap \n";
	param_out << dq.collision_coef      << " \tcollision_coef \n";
	param_out << dq.type_move    << " \ttype_move \n";
	param_out << dq.noise_coef    << " \tnoise_coef \n";
	param_out << dq.use_grid    << " \tuse grid \n";
	//param_out << dq.grid_dim    << " \tgrid dimension \n";
	param_out << dq.dens_cell_size    << " \tcell size for density distribution \n";
	param_out << dq.do_harmonic_forces    << " \tdo harmonic forces \n";
	param_out << dq.harmonic_cnst_x    << " \tharmonic force constant in X\n";
	param_out << dq.harmonic_cnst_y    << " \tharmonic force constant in Y\n";
	param_out << dq.harmonic_cnst_z    << " \tharmonic force constant in Z\n";
	param_out << dq.do_kill    << " \tkill neurons\n";
	param_out << dq.type_kill    << " \ttype_kill\n";
	param_out << dq.kill_prob    << " \tkilling probability\n";
	param_out << dq.kill_bomb_prob    << " \tkilling bomb probability\n";
	param_out << dq.bomb_radius    << " \tkilling bomb radius\n";

	// =============================================
	// FIND BOX SIZE TO INCLUDE TISSUE THICKNESS.
	// THIS IS CORRECTED BELOW, AFTER PUTTING THE NEURONAL
	// COLUMNS.
	// NOTE THAT A SQUARE BOX IS NOT APPROPRIATE FOR PERIODIC BOUNDARY
	// CONDITIONS USING A HEXAGONAL LATTICE. HAVE TO USE RECTANGLE. BUT
	// AS LONG AS THE SLAB IS AWAY FROM THE BOUNDARY THAT IS OK.
	float halfw = dq.tissue_thick / 2.0;
	float halft = dq.L_tissue / 2.0;
	// 		SMALLEST *CIRCLE* TO ENCOMPASS ORIENTATION OF THE BLOCK
	//dq.L_box = 2.0 * sqrt( CUAD(halfw) + CUAD(halft) );
	dq.L_box = dq.L_tissue;
	cout << "1. Current size of system box: " << dq.L_box << "\n";
	//   ACCOMODATE A SLAB WITHOUT CLIPPING OF THE CORNERS WHEN THE CUBE
	//   IS ROTATED TO, e.g., phi=54.833 theta=51.342 
	//dq.L_box /= cos(PI_4);
	dq.L_box = 2.0 * sqrt( CUAD(halfw) + 2.0 * CUAD(halft) );
	cout << "2. Corrected L_box for rotated cube: " << dq.L_box << "\n";

	// ==========================================================
	// GENERATE 3D NEURONAL LANDSCAPE
	float inter_y = 2.0 * sin( RAD(60.) ) * dq.inter_column;

	int num_eveneu_x = int( dq.L_box / dq.inter_column ) + 1;
	int num_eveneu_y = int( dq.L_box / inter_y ) + 1;
	int num_oddneu_x = int( (dq.L_box-(dq.inter_column/2.0)) / dq.inter_column ) + 1;
	int num_oddneu_y = int( (dq.L_box-(dq.inter_column*sin(RAD(60.)))) / inter_y ) + 1;

	// dq.inter_columns STILL DETERMINES TOTAL NO. COLUMNS EVEN IN RANDOM (NOT HEXAGONAL) CASE
	int num_total_columns = num_eveneu_x * num_eveneu_y + num_oddneu_x * num_oddneu_y;
	param_out << num_total_columns << " \ttotal no. columns in ensemble   \n";
	cout << "-- Total no. columns: " << num_total_columns << "\n";

	// MATRICES FOR BASE OF COLUMNS
	float **basemat_even_x;
	float **basemat_even_y;
	float **basemat_odd_x;
	float **basemat_odd_y;
	// ALLOCATE MEMORY
	basemat_even_x = new float*[num_eveneu_x];
	basemat_even_y = new float*[num_eveneu_x];
	for (i=0; i<num_eveneu_x; i++)
	{
	  basemat_even_x[i] = new float[num_eveneu_y];
	  basemat_even_y[i] = new float[num_eveneu_y];
	}
	basemat_odd_x = new float*[num_oddneu_x];
	basemat_odd_y = new float*[num_oddneu_x];
	for (i=0; i<num_oddneu_x; i++)
	{
	  basemat_odd_x[i] = new float[num_oddneu_y];
	  basemat_odd_y[i] = new float[num_oddneu_y];
	}

	// INITIALIZE THE MATRICES FOR CRYSTALLINE PLACEMENT
	float max_x = 0.;
	float max_y = 0.;
	for (i=0; i<num_eveneu_x; i++)
	{
	  for (j=0; j<num_eveneu_y; j++)
	  {
		basemat_even_x[i][j] = i * dq.inter_column;
		basemat_even_y[i][j] = j * inter_y;
		max_x = MAX(max_x, basemat_even_x[i][j] );
		max_y = MAX(max_y, basemat_even_y[i][j] );
	  }
	}
	for (i=0; i<num_oddneu_x; i++)
	{
	  for (j=0; j<num_oddneu_y; j++)
	  {
		basemat_odd_x[i][j] = i * dq.inter_column + dq.inter_column / 2.0;
		basemat_odd_y[i][j] = j * inter_y + dq.inter_column * sin( RAD(60.) );
		max_x = MAX(max_x, basemat_odd_x[i][j] );
		max_y = MAX(max_y, basemat_odd_y[i][j] );
	  }
	}

	// NOW CORRECT SYSTEM SIZE SO THAT COLUMNS AT THE BOUNDARY ARE
	// NOT TOO CLOSE TOGETHER.
	cout << "3. Got for columns: Max x = " << max_x << " Max y = " << max_y << "\n";
	float L_box_tmp = MAX(max_x, max_y) + 2.1*dq.neu_radius;
	dq.L_box = MAX(dq.L_box, L_box_tmp);
	cout << "4. Corrected L_box to separate columns: " << dq.L_box << "\n";
	// WRITE TO PARAM FILE.
	param_out << dq.L_box    << " \tL_box   \n";

	// REASSIGN (OR NOT) COLUMN POSITIONS TO HAVE RANDOM X,Y POSITIONS.
	// RANDOMIZE WITHIN max_x AND max_y.
	char col_pos_pre[40] = "column_xy_";
	strcat(col_pos_pre, file_lab);
	strcat(col_pos_pre, ".xy");
	ofstream col_pos(col_pos_pre);
	if (dq.put_random_col_xy) 
	{
	  cout << "NOTE: putting columns at random positions!\n";

	  // THESE ARE SO THAT COLUMNS ARE ALWAYS ONE RADIUS AWAY FROM THE BORDERS,
	  // RANDOM POSITIONS ARE BETWEEN L_o AND L_o + L_rnd.
	  // COLUMNS WILL AVOID OVERLAP DESCRIBED BY 2*RADIUS_NEURONS.
	  double L_rnd = dq.L_box - 2 * dq.neu_radius;
	  double L_o   = dq.neu_radius;
	  double x1,y1,z1;
	  int ii,jj;
	  int xc, yc;
	  int xcc, ycc;

	  // PUT EVEN
	  for (i=0; i<num_eveneu_x*num_eveneu_y; i++)
	  {
		// GET I,Y INDEX FOR COLUMNS
		yc = int ( i / num_eveneu_x );
		xc = i - yc * num_eveneu_x;

		// GENERATE RANDOM POSITION AND TEST FOR OVERLAP
		int overlap = 1;
		while (overlap)
		{
		  overlap = 0;

		  // GENERATE RANDOM POSITION
		  x1 = L_rnd * rand01();
		  y1 = L_rnd * rand01();
		  z1 = 0.;

		  // SHIFT AWAY FROM WALLS
		  x1+= L_o;
		  y1+= L_o;

		  vec pos(x1, y1, z1);
		  vec diff;

		  // TEST WITH OTHER PLACED COLUMNS
		  for (ii=0; ii<i; ii++)
		  {
			// GET I,Y INDEX FOR COLUMNS
			ycc = int ( ii / num_eveneu_x );
			xcc = ii - ycc * num_eveneu_x;

			vec test( basemat_even_x[xcc][ycc],  basemat_even_y[xcc][ycc], 0. );
			diff = pos - test;

			double r2    = diff.norm2();
			//double crit2 = CUAD(2.*dq.neu_radius);
			double crit2 = CUAD(dq.random_col_mindist);

			// CHECK WITHIN RANGE OF INTERACTION
			//if (r2 <= crit2) overlap = 1;
			if (r2 <= crit2) 
			{
			  // EXCLUDE WITH A HARD WALL
			  // overlap = 1;

			  // EXCLUDE WITH A LINEARLY DECAYING PROBABILITY
		      if (rand01() < ABS(r2-crit2)/crit2) 
				overlap = 1;
			}
		  }
		}

		// POSITION ACCEPTED
		basemat_even_x[xc][yc] = x1;
		basemat_even_y[xc][yc] = y1;
		if (dq.write_files)
			col_pos << x1 << " " << y1 << "\n";
		cout << ".";
	  }
	  // PUT ODD
	  for (i=0; i<num_oddneu_x*num_oddneu_y; i++)
	  {
		// GET I,Y INDEX FOR COLUMNS
		yc = int ( i / num_oddneu_x );
		xc = i - yc * num_oddneu_x;

		// GENERATE RANDOM POSITION AND TEST FOR OVERLAP
		int overlap = 1;
		while (overlap)
		{
		  overlap = 0;

		  // GENERATE RANDOM POSITION
		  x1 = L_rnd * rand01();
		  y1 = L_rnd * rand01();
		  z1 = 0.;

		  // SHIFT AWAY FROM WALLS
		  x1+= L_o;
		  y1+= L_o;

		  vec pos(x1, y1, z1);
		  vec diff;

		  // TEST WITH OTHER PLACED COLUMNS
		  // -- EVEN FIRST (ALL OF THEM)
		  for (ii=0; ii<num_eveneu_x*num_eveneu_y; ii++)
		  {
			// GET I,Y INDEX FOR COLUMNS
			ycc = int ( ii / num_eveneu_x );
			xcc = ii - ycc * num_eveneu_x;

			vec test( basemat_even_x[xcc][ycc],  basemat_even_y[xcc][ycc], 0. );
			diff = pos - test;

			double r2    = diff.norm2();
			//double crit2 = CUAD(2.*dq.neu_radius);
			double crit2 = CUAD(dq.random_col_mindist);

			// CHECK WITHIN RANGE OF INTERACTION
			if (r2 <= crit2) 
			{
			  // EXCLUDE WITH A HARD WALL
			  // overlap = 1;

			  // EXCLUDE WITH A LINEARLY DECAYING PROBABILITY
		      if (rand01() < ABS(r2-crit2)/crit2) 
				overlap = 1;
			}
		  }
		  // -- ODD NOW 
		  for (ii=0; ii<i; ii++)
		  {
			// GET I,Y INDEX FOR COLUMNS
			ycc = int ( ii / num_oddneu_x );
			xcc = ii - ycc * num_oddneu_x;

			vec test( basemat_odd_x[xcc][ycc],  basemat_odd_y[xcc][ycc], 0. );
			diff = pos - test;

			double r2    = diff.norm2();
			//double crit2 = CUAD(2.*dq.neu_radius);
			double crit2 = CUAD(dq.random_col_mindist);

			// CHECK WITHIN RANGE OF INTERACTION
			if (r2 <= crit2) 
			{
			  // EXCLUDE WITH A HARD WALL
			  // overlap = 1;

			  // EXCLUDE WITH A LINEARLY DECAYING PROBABILITY
		      if (rand01() < ABS(r2-crit2)/crit2) 
				overlap = 1;
			}
		  }
		}

		// POSITION ACCEPTED
		basemat_odd_x[xc][yc] = x1;
		basemat_odd_y[xc][yc] = y1;
		if (dq.write_files)
			col_pos << x1 << " " << y1 << "\n";
		cout << ".";
	  }
	  cout << "\n";
	  cout << "DONE: putting columns at random positions!\n";
	}

	// INITIALIZE GRID
	grid neu_grid(dq.L_box, dq.grid_dim, dq.neu_radius);
	param_out << neu_grid.cell_size    << " \tsize of grid   \n";
	cout << "2*radius = " << 2.*dq.neu_radius << "\n";

	// CONSTRUCT THE WHOLE NEURONAL CUBIC LANDSCAPE
	infobjarray<neuron,NEU_ARRAY> neural_array; 
	float z_neu_init, z_neu_init_0;
	int   num_neu_init = 0;

	// EXAMPLE COLUMN TO BE MARKED FOR WRITING
	//int marked_col_x_id = int(num_eveneu_x / 2.);
	//int marked_col_y_id = int(num_eveneu_y / 2.);
	int marked_col_to_write = int( (num_eveneu_x+1) * num_eveneu_y / 2. );
	// - NOW, ALL NEURONS ARE MARKED BY THE ID (>0) OF THEIR COLUMN -
	// - INTERSTITIAL NEURONS HAVE ID = -1                          -
	int marked_id = 1;

	cout << " [Constructing neuronal landscape.]\n";
	//  DO EVEN
	for (i=0; i<num_eveneu_x; i++)
	{
	  for (j=0; j<num_eveneu_y; j++)
	  {
		// GROW COLUMN UPWARDS
		if (dq.put_nice_z)
		  z_neu_init   = dq.L_tissue/2.;
		else
		  z_neu_init   = dq.L_tissue/2. + random_dev(dq.L_tissue/2.);
		z_neu_init_0 = z_neu_init ;
		// FILL UP TO THE TOP, L_box, WITH NEURONS
		//   ALLOW FOR DEVIATIONS IN THE XY POS OF THE COLUMNS
		double ranx;
		double rany;
		// CHECK FOR CRYSTALLINE OR RANDOM COLUMN PLACEMENT
		if (!dq.put_random_col_xy) 
		{
		  ranx = basemat_even_x[i][j] + random_dev(dq.dev_base);
		  rany = basemat_even_y[i][j] + random_dev(dq.dev_base);
		  if (dq.write_files)
		  	col_pos << ranx << " " << rany << "\n";
		}
		else
		{
		  ranx = basemat_even_x[i][j];
		  rany = basemat_even_y[i][j];
		  // THIS ONE IS WRITEN TO FILE SOMEWHERE ELSE
		}

		while (z_neu_init < dq.L_box)
		{
		  // SKIP skip_neurons PERCENT OF NEURONS
		  if (rand01()>dq.skip_neurons)
		  {
			neuron *temp_neuron = new neuron;
			temp_neuron->setOriginalPosition( vec( ranx + random_dev(dq.delta_xy_neuron) , rany + random_dev(dq.delta_xy_neuron) , z_neu_init ) );
			temp_neuron->setRadius(dq.neu_radius);
			temp_neuron->setEven();

			// IF CHOSEN FOR MARKING
			//if ( i==marked_col_x_id && j==marked_col_y_id ) temp_neuron->setMarked();
			// MARK WITH ID OF COLUMN
			temp_neuron->setMarked(marked_id);

			neural_array.add_one(*temp_neuron);
			num_neu_init++;
		  }
		  //z_neu_init += dq.inter_neuron + random_dev(dq.delta_neuron);
		  z_neu_init += dq.inter_neuron + gaussian_dist() * dq.delta_neuron;
		}
		// GROW COLUMN DOWNWARDS
		z_neu_init  = z_neu_init_0 ;
		if (dq.put_nice_z)
		  z_neu_init -= dq.inter_neuron;
		else
		  z_neu_init -= dq.inter_neuron + random_dev(dq.delta_neuron);
		while (z_neu_init > 0)
		{
		  // SKIP skip_neurons PERCENT OF NEURONS
		  if (rand01()>dq.skip_neurons)
		  {
			neuron *temp_neuron = new neuron;
			temp_neuron->setOriginalPosition( vec( ranx + random_dev(dq.delta_xy_neuron) , rany + random_dev(dq.delta_xy_neuron) , z_neu_init ) );
			temp_neuron->setRadius(dq.neu_radius);
			temp_neuron->setEven();

			// IF CHOSEN FOR MARKING
			//if ( i==marked_col_x_id && j==marked_col_y_id ) temp_neuron->setMarked();
			// MARK WITH ID OF COLUMN
			temp_neuron->setMarked(marked_id);

			neural_array.add_one(*temp_neuron);
			num_neu_init++;
		  }
		  //z_neu_init -= dq.inter_neuron + random_dev(dq.delta_neuron);
		  z_neu_init -= dq.inter_neuron + gaussian_dist() * dq.delta_neuron;
		}

		// INCREASE ID OF COLUMN
		marked_id++;

		// DONE WITH THIS COLUMN
	  }
	}

	//  DO ODD
	for (i=0; i<num_oddneu_x; i++)
	{
	  for (j=0; j<num_oddneu_y; j++)
	  {
		// GROW COLUMN UPWARDS
		if (dq.put_nice_z)
		  z_neu_init   = dq.L_tissue/2.;
		else
		  z_neu_init   = dq.L_tissue/2. + random_dev(dq.L_tissue/2.);
		z_neu_init_0 = z_neu_init ;
		// FILL UP TO THE TOP, L_box, WITH NEURONS
		//   ALLOW FOR DEVIATIONS IN THE XY POS OF THE COLUMNS
		double ranx;
		double rany;
		// CHECK FOR CRYSTALLINE OR RANDOM COLUMN PLACEMENT
		if (!dq.put_random_col_xy) 
		{
		  ranx = basemat_odd_x[i][j] + random_dev(dq.dev_base);
		  rany = basemat_odd_y[i][j] + random_dev(dq.dev_base);
		  if (dq.write_files)
		  	col_pos << ranx << " " << rany << "\n";
		}
		else
		{
		  ranx = basemat_odd_x[i][j];
		  rany = basemat_odd_y[i][j];
		  // THIS ONE IS WRITEN TO FILE SOMEWHERE ELSE
		}

		while (z_neu_init < dq.L_box)
		{
		  // SKIP skip_neurons PERCENT OF NEURONS
		  if (rand01()>dq.skip_neurons)
		  {
			neuron *temp_neuron = new neuron;
			temp_neuron->setOriginalPosition( vec( ranx + random_dev(dq.delta_xy_neuron), rany + random_dev(dq.delta_xy_neuron) , z_neu_init ) );
			temp_neuron->setRadius(dq.neu_radius);
			temp_neuron->setOdd();

			// MARK WITH ID OF COLUMN
			temp_neuron->setMarked(marked_id);

			neural_array.add_one(*temp_neuron);
			num_neu_init++;
		  }
		  //z_neu_init += dq.inter_neuron + random_dev(dq.delta_neuron);
		  z_neu_init += dq.inter_neuron + gaussian_dist() * dq.delta_neuron;
		}
		// GROW COLUMN DOWNWARDS
		z_neu_init   = z_neu_init_0 ;
		if (dq.put_nice_z)
		  z_neu_init -= dq.inter_neuron;
		else
		  z_neu_init -= dq.inter_neuron + random_dev(dq.delta_neuron);
		while (z_neu_init > 0)
		{
		  // SKIP skip_neurons PERCENT OF NEURONS
		  if (rand01()>dq.skip_neurons)
		  {
			neuron *temp_neuron = new neuron;
			temp_neuron->setOriginalPosition( vec( ranx + random_dev(dq.delta_xy_neuron), rany + random_dev(dq.delta_xy_neuron) , z_neu_init ) );
			temp_neuron->setRadius(dq.neu_radius);
			temp_neuron->setOdd();

			// MARK WITH ID OF COLUMN
			temp_neuron->setMarked(marked_id);

			neural_array.add_one(*temp_neuron);
			num_neu_init++;
		  }
		  //z_neu_init -= dq.inter_neuron + random_dev(dq.delta_neuron);
		  z_neu_init -= dq.inter_neuron + gaussian_dist() * dq.delta_neuron;
		}

		// INCREASE ID OF COLUMN
		marked_id++;

		// DONE WITH THIS COLUMN
	  }
	}

	// CALCULATE PAIR CORRELATION FUNCTION.
	pair_corr(num_eveneu_x, num_eveneu_y, num_oddneu_x, num_oddneu_y, basemat_even_x, basemat_even_y, basemat_odd_x, basemat_odd_y, file_lab, dq);

	// CLOSE FILE OF COLUMN XY POSITIONS
	col_pos.close();

	// PUT INTERSTICIAL NEURONS RANDOMLY WHEREVER THEY FIT.
	cout << " [Putting intersticial neurons in empty cavities.]\n";
	//  * THE NUMBER OF INTERSTICIAL NEURONS IS dq.interst_neurons PERCENT
	//  * OF THE *TOTAL* NUMBER OF NEURONS. THEN THE TOTAL NUMBER OF THESE
	//  * IS GIVEN BY:    n * P / (1-P)
	int tot_interst = ANINT( (num_neu_init * dq.interst_neurons) / (1.0 - dq.interst_neurons) );
	cout << "   Will put " << tot_interst << " intersticial neurons, \n"
	  << "   on top of " << num_neu_init << " lattice neurons,\n"
	  << "   for a total of " << num_neu_init + tot_interst 
	  << "   neurons.\n";
	int how_many_interst = put_interst_neurons(neural_array, dq, tot_interst);
	num_neu_init += how_many_interst;
	float real_interst = how_many_interst*100./tot_interst;
	cout << "   Placed " << real_interst << "% (" << how_many_interst << ") "
	  << "of the intersticial neurons\n";
	param_out << dq.interst_neurons     << " \tintersticial neurons (fraction)\n";
	param_out << how_many_interst*1.0/num_neu_init << " \tActual intersticial (fraction)   \n";
	param_out << how_many_interst  << " \tActual intersticial no.   \n";

	// SET TOTAL NUMBER OF PARTICLES.
	dq.num_neu = num_neu_init;
	cout << "Generated " << dq.num_neu << " initial neurons. \n";
	param_out << dq.num_neu    << " \tinitial neurons   \n";
	param_out << dq.num_neu/(dq.L_box*dq.L_box*dq.L_box) << " \tinitial neu/pix^3   \n";

	// ==========================================================
	// KILL NEURONS ACCORDING TO TYPE OF KILLING.
	if (dq.do_kill)
	{
	  cout << " [Killing neurons.]\n";

	  int how_many_killed = kill_neurons(neural_array, dq, file_lab);
	  // UPDATE TOTAL NUM OF NEURONS
	  dq.num_neu -= how_many_killed;
	  // GET CURRENT TOTAL NUM OF INTERSTICIAL NEURONS
	  int interst = interst_neurons(neural_array);
	  double interst_per = interst*1.0/dq.num_neu;

	  cout << "Neurons after killing: " << dq.num_neu << "\n";
	  cout << "Interst after killing: " << interst << "\n";
	  cout << "Actual Interst % after killing: " << interst_per * 100. << "\n";

	  param_out << dq.num_neu    << " \tafterkill total neurons   \n";
	  param_out << dq.num_neu/(dq.L_box*dq.L_box*dq.L_box) << " \tafterkill neu/pix^3   \n";
	  param_out << interst_per    << " \tafterkill Actual intersticials (fraction)  \n";
	  param_out << interst        << " \tafterkill Actual intersticials no.\n";
	} 

	// ==========================================================
	// RELAX THE SYSTEM ELIMINATING ALL INITIAL COLLISIONS AND OVERLAPS
	dq.do_relax = 1;
	cout << " [Relaxing system to avoid overlap...]\n";
	//  * WRITE BLOCK BEFORE RELAXING.
	if (dq.write_files)
	{
	  char pre_relax[40]     = "preRelaxBlock_01_";
	  strcat(pre_relax, file_lab);
	  strcat(pre_relax, ".xyz");
	  write_block_to_file(neural_array, dq, pre_relax);
	}
	//  * RELAX
	int num_relax = 0;
	do {
	  //cout << "to repopulate\n";
	  neu_grid.repopulate(neural_array,dq);
	  //cout << "to calcandmove\n";
	  calc_and_move(neural_array, neu_grid, dq, 0, col);
	  //cout << "to *\n";
	  num_relax++;
	  cout << "t_relax = " << num_relax << "\n";
	} while( !neuralForcesZero(neural_array) );
	cout << "Relaxed system by " << num_relax << " steps.\n";
	param_out << num_relax    << " \tRelaxing steps \n";
	//  * UPDATE ORIGINAL POSITIONS
	updateOriginalPositions(neural_array);

	// WRITE THE 3D BLOCK OF PRISTINE BRAIN TO FILE.
	//  * ONLY WRITE IF NOT WRITING FILES vs TIME.
	//  * THIS INIT FILE IS THE SAME AS THE t=0 FILE.
	if (dq.write_files)
	{
	  char pre_file[40]     = "preMoveBlock_02_";
	  strcat(pre_file, file_lab);
	  strcat(pre_file, ".xyz");
	  //
	  write_block_to_file(neural_array, dq, pre_file);

	  // WRITE INTERSTICIAL NEURONS
	  char pre_interst[40]     = "interst_";
	  strcat(pre_interst, file_lab);
	  strcat(pre_interst, ".xyz");
	  //
	  ofstream interst_neu(pre_interst);
	  for (i=0; i<=neural_array.last_one; i++)
	  {
		// WRITE TO FILE THE WRAPPED POSITION INSIDE THE BOX.
		if ( neural_array[i].isAlive() && neural_array[i].isInterst())
		  interst_neu << neural_array[i].getWrapPosition(dq.L_box) << "\n";
	  }
	  interst_neu.close();

	  // WRITE DEAD NEURONS
	  char dead[40]     = "dead_";
	  strcat(dead, file_lab);
	  strcat(dead, ".xyz");
	  //
	  ofstream dead_neu(dead);
	  for (i=0; i<=neural_array.last_one; i++)
	  {
		// WRITE TO FILE THE WRAPPED POSITION INSIDE THE BOX.
		if ( !neural_array[i].isAlive() )
		  dead_neu << neural_array[i].getWrapPosition(dq.L_box) << "\n";
	  }
	  dead_neu.close();

	  // WRITE MATRIX OF NEURONAL DENSITIES
	  char dens_map[40]     = "densMapCells_";
	  strcat(dens_map, file_lab);
	  strcat(dens_map, ".dat");
	  //
	  write_density_distribution(neural_array, dq, dens_map);

	  // WRITE MARKED COLUMN
	  char marked[40] = "marked_column_";
	  strcat(marked, file_lab);
	  strcat(marked, ".xyz");
	  ofstream marked_neu(marked);
	  for (i=0; i<=neural_array.last_one; i++)
	  {
		// WRITE TO FILE THE WRAPPED POSITION INSIDE THE BOX.
		if ( neural_array[i].getMarked() == marked_col_to_write )
		  marked_neu << neural_array[i].getWrapPosition(dq.L_box) << "\n";
	  }
	  marked_neu.close();

	  // WRITE BACKGROUND FOR THE MARKED COLUMN
	  char bck_sys[40] = "background_system_";
	  strcat(bck_sys, file_lab);
	  strcat(bck_sys, ".xyz");
	  ofstream bck_sys_neu(bck_sys);
	  // * FIRST, GET A MARKED NEURON
	  vec which;
	  for (i=0; i<=neural_array.last_one; i++)
	  {
		if ( neural_array[i].getMarked() == marked_col_to_write )
		  which = neural_array[i].getWrapPosition(dq.L_box);
	  }
	  cout << "got marked neuron " << which << "\n";
	  double xmn = which.x;
	  double ymn = which.y;
	  // * NOW WRITE ALL NEURONS EXCEPT THOSE IN COLUMNS IN FRONT OF MARKED
	  for (i=0; i<=neural_array.last_one; i++)
	  {
		if ( neural_array[i].isAlive() )
		{
		  if ( neural_array[i].getWrapPosition(dq.L_box).x < xmn )
		  {
			bck_sys_neu << neural_array[i].getWrapPosition(dq.L_box) << "\n";
		  }
		}
	  }
	  bck_sys_neu.close();

	}


	// ==========================================================
	// MOVE THE NEURONS BY THE SPECIFIED NUMBER OF STEPS
	cout << " [Moving neurons.]\n";

	//  * NOT A RELAXING STEP.
	dq.do_relax = 0;
	//  * RESET COLLISION ARRAY
	col.reset();

	char pre_r[40]     = "R_meandisp_";
	strcat(pre_r, file_lab);
	strcat(pre_r, ".dat");
	ofstream r_mean(pre_r);
	char pre_r_XY[40]     = "R_XY_meandisp_";
	strcat(pre_r_XY, file_lab);
	strcat(pre_r_XY, ".dat");
	ofstream r_mean_XY(pre_r_XY);
	char pre_r_Z[40]     = "R_Z_meandisp_";
	strcat(pre_r_Z, file_lab);
	strcat(pre_r_Z, ".dat");
	ofstream r_mean_Z(pre_r_Z);

	double ave_Rdis, stdev_Rdis;
	double ave_Rdis_XY, stdev_Rdis_XY;
	double ave_Rdis_Z, stdev_Rdis_Z;
	for (i=0; i<=dq.num_steps; i++)
	{
	  // MEAN DISPLACEMENT
	  mean_neural_displacement(
		  ave_Rdis, stdev_Rdis, 
		  ave_Rdis_XY, stdev_Rdis_XY, 
		  ave_Rdis_Z, stdev_Rdis_Z, 
		  neural_array, dq, dq.L_box
		  );
	  cout << "t = " << i
		<< " : Average displacement is " 
		<< ave_Rdis  << " +- " << stdev_Rdis
		<< " : XY displacement is " 
		<< ave_Rdis_XY  << " +- " << stdev_Rdis_XY
		<< " : Z displacement is " 
		<< ave_Rdis_Z  << " +- " << stdev_Rdis_Z
		<< "\n";
	  if (dq.write_files)
	  {
		r_mean << i << " " << ave_Rdis << " " << stdev_Rdis << "\n";
		r_mean_XY << i << " " << ave_Rdis_XY << " " << stdev_Rdis_XY << "\n";
		r_mean_Z << i << " " << ave_Rdis_Z << " " << stdev_Rdis_Z << "\n";
	  }

	  // WRITE NEURONAL BLOCK PERIODICALLY TO FILE
	  if ( ((i % dq.write_freq ) == 0) && dq.write_periodic_files )
	  {
		// WRITE THE MOVED 3D BLOCK
		char pre_f[40]     = "neuBlock_";
		strcat(pre_f, file_lab);
		strcat(pre_f, "_");
		char *cpix_file     = NULL;
		cpix_file     = conc_char_2_int(pre_f,i,3);
		strcat(cpix_file,     ".xyz");

		write_block_to_file(neural_array, dq, cpix_file);
		delete [] cpix_file;
	  }

	  // MOVE	
	  if (i<dq.num_steps)
	  {
#if DEBUG == 1
		cout << " * to repo " << i << "\n";
#endif

		neu_grid.repopulate(neural_array,dq,i);

#if DEBUG == 1
		cout << " * to move " << i << "\n";
#endif

		calc_and_move(neural_array, neu_grid, dq, i, col);

#if DEBUG == 1
		cout << " * out of move " << i << "\n";
#endif
	  }
	}
	  if (dq.write_files)
	  {
		r_mean.close();
		r_mean_XY.close();
		r_mean_Z.close();
	  }

	// CALCULATE MEAN DISPLACEMENT DISTANCE
	mean_neural_displacement(
		ave_Rdis, stdev_Rdis, 
		ave_Rdis_XY, stdev_Rdis_XY, 
		ave_Rdis_Z, stdev_Rdis_Z, 
		neural_array, dq, dq.L_box
		);
	//cout << "Average displacement is " << ave_Rdis << "\n";
	param_out << ave_Rdis    << " \tAverage displacement   \n";
	param_out << stdev_Rdis  << " \tStd. Dev. displacement   \n";
	param_out << ave_Rdis_XY    << " \tAverage XY displacement   \n";
	param_out << stdev_Rdis_XY  << " \tStd. Dev. XY displacement   \n";
	param_out << ave_Rdis_Z    << " \tAverage Z displacement   \n";
	param_out << stdev_Rdis_Z  << " \tStd. Dev. Z displacement   \n";

	// WRITE DISTRIBUTION OF DISPLACEMENTS
	if (dq.write_files)
	{
	  char disp_dist[40]     = "dist_distri_";
	  strcat(disp_dist, file_lab);
	  strcat(disp_dist, ".dat");
	  write_displacement_distribution(neural_array, disp_dist);

	  // WRITE COLLISIONS DURING MOVING
	  char col_file[40]     = "collisions_";
	  strcat(col_file, file_lab);
	  strcat(col_file, ".dat");
	  col.write(col_file);

	  // WRITE A GIVEN CELL
	  //int a=0,b=1,c=1;
	  //neu_grid.write_cell(a,b,c,file_lab,neural_array,dq);
	}

	// WRITE FINAL BLOCK OF MOVED NEURONS WITH COLUMN ID
	if (dq.write_files)
	{
		// WRITE THE FINALIZED MOVED 3D BLOCK
		char pre_f[40]     = "neuBlockwithColID_";
		strcat(pre_f, file_lab);
		strcat(pre_f, ".xyz");

		write_block_to_file_colID(neural_array, dq, pre_f);
	}

	// ==========================================================
	// PUT GLIA INSIDE. THIS ASSUMES THAT GLIA DOES NOT INTERFERE
	// WITH NEURONAL DISPLACEMENTS, THUS CAN BE PLACED AT THE END.
	//
	// CONSTRUCT THE GLIA ARRAY
	infobjarray<glia,GLIA_ARRAY> glial_array; 

	if (dq.put_glia == 1)
	{
	  cout << " [Putting glia in empty spaces.]\n";
	  int tot_glia = dq.glia_amount * dq.num_neu;
	  cout << "   Will put " << tot_glia << " glia. \n";
	  int how_many_glia = put_glia(glial_array, neural_array, dq, tot_glia);
	  float real_glia = how_many_glia*100./tot_glia;	
	  cout << "   Placed " << real_glia << "% (" << how_many_glia << ") "
		<< "of the glia\n";
	  param_out << dq.glia_amount     << " \tglia (fraction)\n";
	  param_out << how_many_glia*1.0/dq.num_neu << " \tActual glia (fraction)   \n";
	  param_out << how_many_glia  << " \tActual glia no.   \n";

	  // SET TOTAL NUMBER OF GLIA.
	  dq.num_glia = how_many_glia;
	  cout << "Generated " << dq.num_glia << " initial glia. \n";
	  param_out << dq.num_glia    << " \tinitial glia   \n";
	  param_out << dq.num_glia/(dq.L_box*dq.L_box*dq.L_box) << " \tinitial glia/pix^3   \n";
	}


	// ==========================================================
	// ROTATE BLOCK OF NEURONS AND GLIA TO ALLOW FOR SIDEWAYS CUT, NOT PARALLEL 
	// TO MICROCOLUMNS
	rotate_neuralblock(neural_array, dq);
	// NOW GLIA
	if (dq.put_glia == 1)
		rotate_glialblock(glial_array, dq);

	// WRITE THE ROTATED 3D NEURONAL BLOCK
	if (dq.write_files)
	{
	  char rot_f[50]     = "postMovedRotNeuralBlock_03_";
	  strcat(rot_f, file_lab);
	  strcat(rot_f, ".xyz");

	  ofstream rot_neu(rot_f);
	  for (j=0; j<=neural_array.last_one; j++)
	  {
		if (neural_array[j].isAlive())
		  rot_neu << neural_array[j].getRotatedPosition() << "\n";
	  }
	  rot_neu.close();
	}

	// WRITE THE ROTATED 3D GLIAL BLOCK
	if (dq.write_files && dq.put_glia)
	{
	  char rot_f[50]     = "postMovedRotGlialBlock_03_";
	  strcat(rot_f, file_lab);
	  strcat(rot_f, ".xyz");

	  ofstream rot_glia(rot_f);
	  for (j=0; j<=glial_array.last_one; j++)
	  {
		if (glial_array[j].isAlive())
		  rot_glia << glial_array[j].getRotatedPosition() << "\n";
	  }
	  rot_glia.close();
	}
	// END: ROTATE BLOCK 

	// ==========================================================
	// 1. FIND NEURONS INSIDE THE TISSUE SLAB, DEFINED BY AN x-y POLYGON
	cout << " [Finding neurons inside slab.]\n";
	// 1.1 GENERATE POINTS FOR THE POLYGON THAT DEFINE THE TISSUE
	float *x_cp = new float[n+1];
	float *y_cp = new float[n+1];
	//  * NOTE: THIS dq.angle IS THE ANGLE OF THE CUTTING KNIFE, NOT THE BLOCK ROTATING ANGLES.
	x_cp[0] = dq.L_box/2.0 + halfw * cos(PI_2 - dq.angle) + halft * cos(dq.angle);
	y_cp[0] = dq.L_box/2.0 - halfw * sin(PI_2 - dq.angle) + halft * sin(dq.angle);
	x_cp[1] = dq.L_box/2.0 - halfw * cos(PI_2 - dq.angle) + halft * cos(dq.angle);
	y_cp[1] = dq.L_box/2.0 + halfw * sin(PI_2 - dq.angle) + halft * sin(dq.angle);
	x_cp[2] = dq.L_box/2.0 - halfw * cos(PI_2 - dq.angle) - halft * cos(dq.angle);
	y_cp[2] = dq.L_box/2.0 + halfw * sin(PI_2 - dq.angle) - halft * sin(dq.angle);
	x_cp[3] = dq.L_box/2.0 + halfw * cos(PI_2 - dq.angle) - halft * cos(dq.angle);
	y_cp[3] = dq.L_box/2.0 - halfw * sin(PI_2 - dq.angle) - halft * sin(dq.angle);
	// COPY FOR CYCLIC CONVERSION.
	x_cp[n] = x_cp[0];
	y_cp[n] = y_cp[0];

	// 1.2 FIND NEURONS OUT vs NEURONS IN
	data_array x_intercpts;
	x_intercpts.init(n);

	// 		LOOP OVER NEURONS AND DETERMINE IF INSIDE THE POLYGON.
	char size_d[40] = "size_dist_";
	strcat(size_d, file_lab);
	strcat(size_d, ".dat");
	ofstream size_dist(size_d);
	//
	int num_interior_neurons = 0;
	for (j=0; j<=neural_array.last_one; j++)
	{
	  if (neural_array[j].isAlive())
	  {
		// LOOP OVER PAIRS OF Y VALUES
		for (k=0; k<n; k++)
		{
		  if ( in_between(neural_array[j].getRotatedPosition().y, y_cp[k], y_cp[k+1]) )
		  {
			float x_0 = x_cp[k+1] - ( x_cp[k+1] - x_cp[k] ) * ( y_cp[k+1] - neural_array[j].getRotatedPosition().y ) / ( y_cp[k+1] - y_cp[k] );
			x_intercpts.insert_value(x_0);
		  }
		}

		// CHECK THAT THE X VALUE IS ALSO IN
		int positives = 0;
		for (k=0; k<x_intercpts.curr_index; k++)
		{
		  if ( (x_intercpts.value[k] - neural_array[j].getRotatedPosition().x) > 0 ) positives++;
		}

		// CHECK Z IS WITHIN LIMITS
		if ( (positives%2)!=0 && (neural_array[j].getRotatedPosition().z < dq.L_tissue) && (neural_array[j].getRotatedPosition().z >= 0.) ) 
		{
		  // FOUND INSIDE 
		  neural_array[j].setInsideSlab();
		  num_interior_neurons++;

		  // -- SIZE DISTRIBUTION --
		  // WRITE OUT THE MINIMUM DISTANCE FROM THE CENTER TO THE CUTTING SURFACE AT THE Y-Z PLANE.
		  // THIS IS TO CALCULATE DISTRIBUTION OF "AREAS" OF NEURONS IN THE SLICE INCLUDING THOSE 
		  // AREAS THAT GET CUT FROM THE DISSECTION OF THE SLICE.
		  float min_y_surf = dq.L_box;
		  float dist_y_surf;
		  float eff_r;
		  for (k=0; k<n; k++)
		  {
			dist_y_surf = ABS( neural_array[j].getRotatedPosition().y - y_cp[k] );
			min_y_surf = ( min_y_surf < dist_y_surf ? min_y_surf : dist_y_surf);
		  }
		  // WRITE OUT THE EFFECTIVE AREA SEEN THROUGH THE MICROSCOPE.
		  //  ** WARNING - FOR SIMPLICITY THIS DISTRIBUTION IS NOT FOR THE CURRENT THICKNESS
		  //               OF SLICE RUNNING HERE.
		  // THE DISTRIBUTION CALCULATED CORRESPONDS TO THAT SEEN IN A SLIDE OF WIDTH
		  // dq.tissue_thick - 2 * dq.neu_radius
		  // IF WANT THE DISTRIBUTION FOR THE CURRENTLY RUN SLIDE WIDTH, RE-RUN WITH CORRESPONDING
		  // BIGGER THICKNESS.
		  // DIFFERENCE COMES IN THE SOLUTION OF THE CIRCLE, GIVING CROSS SECTION BELOW AND NOT
		  // ABOVE CENTER OF NEURON.
		  if (dq.write_files)
		  {
			if ( min_y_surf > dq.neu_radius )
			{
			  //cout << "area " << 4.*PI*CUAD(dq.neu_radius) << "\n";
			  size_dist << 4.*PI*CUAD(dq.neu_radius) << "\n";
			}
			else
			{
			  eff_r = sqrt( CUAD(dq.neu_radius) - CUAD(dq.neu_radius - min_y_surf) );
			  //cout << "area " << 4.*PI*CUAD(eff_r) << "\n";
			  size_dist << 4.*PI*CUAD(eff_r) << "\n";
			}
		  }
		  //cout << "neuron " << min_y_surf << " away from y-z plane\n";

		}
		x_intercpts.init(n);

	  }
	}
	size_dist.close();
	cout << "Found " << num_interior_neurons << " neurons in the slab.\n";
	param_out << num_interior_neurons    << " \tneurons in the slab   \n";

	// 2. NOW, FIND GLIA INSIDE THE TISSUE SLAB, DEFINED BY AN x-y POLYGON
	if (dq.put_glia == 1)
	{
	  cout << " [Finding glia inside slab.]\n";

	  // 2.1 FIND GLIA OUT vs GLIA IN
	  x_intercpts.init(n);

	  // 		LOOP OVER GLIA AND DETERMINE IF INSIDE THE POLYGON.
	  int num_interior_glia = 0;
	  for (j=0; j<=glial_array.last_one; j++)
	  {
		if (glial_array[j].isAlive())
		{
		  // LOOP OVER PAIRS OF Y VALUES
		  for (k=0; k<n; k++)
		  {
			if ( in_between(glial_array[j].getRotatedPosition().y, y_cp[k], y_cp[k+1]) )
			{
			  float x_0 = x_cp[k+1] - ( x_cp[k+1] - x_cp[k] ) * ( y_cp[k+1] - glial_array[j].getRotatedPosition().y ) / ( y_cp[k+1] - y_cp[k] );
			  x_intercpts.insert_value(x_0);
			}
		  }

		  // CHECK THAT THE X VALUE IS ALSO IN
		  int positives = 0;
		  for (k=0; k<x_intercpts.curr_index; k++)
		  {
			if ( (x_intercpts.value[k] - glial_array[j].getRotatedPosition().x) > 0 ) positives++;
		  }

		  // CHECK Z IS WITHIN LIMITS
		  if ( (positives%2)!=0 && (glial_array[j].getRotatedPosition().z < dq.L_tissue) && (glial_array[j].getRotatedPosition().z >= 0.) ) 
		  {
			// FOUND INSIDE 
			glial_array[j].setInsideSlab();
			num_interior_glia++;
		  }
		  x_intercpts.init(n);
		}
	  }
	  cout << "Found " << num_interior_glia << " glia in the slab.\n";
	  param_out << num_interior_glia    << " \tglia in the slab   \n";
	}
	// END: FIND THINGS INSIDE THE TISSUE SLAB, DEFINED BY AN x-y POLYGON

	// ==========================================================
	//  * * * * * START: ONLY FOR BOMB KILLINGS * * * * * * 
	// 		LOOP OVER BOMB LOCATIONS AND DETERMINE IF INSIDE THE POLYGON.
	x_intercpts.init(n);
	for (j=0; j<=neural_array.last_one; j++)
	{
	  if (neural_array[j].isBombLoc())
	  {
		// LOOP OVER PAIRS OF Y VALUES
		for (k=0; k<n; k++)
		{
		  if ( in_between(neural_array[j].getRotatedPosition().y, y_cp[k], y_cp[k+1]) )
		  {
			float x_0 = x_cp[k+1] - ( x_cp[k+1] - x_cp[k] ) * ( y_cp[k+1] - neural_array[j].getRotatedPosition().y ) / ( y_cp[k+1] - y_cp[k] );
			x_intercpts.insert_value(x_0);
		  }
		}

		// CHECK THAT THE X VALUE IS ALSO IN
		int positives = 0;
		for (k=0; k<x_intercpts.curr_index; k++)
		{
		  if ( (x_intercpts.value[k] - neural_array[j].getRotatedPosition().x) > 0 ) positives++;
		}

		// CHECK Z IS WITHIN LIMITS
		if ( (positives%2)!=0 && (neural_array[j].getRotatedPosition().z < dq.L_tissue) && (neural_array[j].getRotatedPosition().z >= 0.) ) 
		{
		  // FOUND INSIDE 
		  neural_array[j].setInsideSlab();
		}
		x_intercpts.init(n);
	  }
	}
	//  * * * * * END: ONLY FOR BOMB KILLINGS * * * * * * 

	// ==========================================================
	// A. WRITE THE 3D BLOCK NEURONS INSIDE THE SLAB TO FILE.
	// ALSO, TRANSFORM FROM THE X-Y-Z BASIS TO THE X'-Y'-Z' ROTATED BASIS.
	// GENERATE X'-Z' DATA ONLY (TISSUE COLLAPSE) WHERE Z=Z'.
	cout << " [Projecting neurons and writing out.]\n";
	//if (dq.write_files)
	//{
	char pre_file2[40]     = "slab_2D_";
	char pre_file3[40]     = "slab_3D_";
	char pre_file4[40]     = "slab_2D_NotInterst_";
	char pre_file5[40]     = "slab_2D_Interst_";
	char pre_file6[40]     = "slab_3D_BombLoc_";
	strcat(pre_file2, file_lab);
	strcat(pre_file2, "_");
	strcat(pre_file3, file_lab);
	strcat(pre_file3, "_");
	strcat(pre_file4, file_lab);
	strcat(pre_file4, "_");
	strcat(pre_file5, file_lab);
	strcat(pre_file5, "_");
	strcat(pre_file6, file_lab);
	strcat(pre_file6, "_");
	char *cpixel_file2     = NULL;
	cpixel_file2     = conc_char_2_int(pre_file2,dq.num_steps,3);
	char *cpixel_file3     = NULL;
	cpixel_file3     = conc_char_2_int(pre_file3,dq.num_steps,3);
	char *cpixel_file4     = NULL;
	cpixel_file4     = conc_char_2_int(pre_file4,dq.num_steps,3);
	char *cpixel_file5     = NULL;
	cpixel_file5     = conc_char_2_int(pre_file5,dq.num_steps,3);
	char *cpixel_file6     = NULL;
	cpixel_file6     = conc_char_2_int(pre_file6,dq.num_steps,3);
	strcat(cpixel_file2,     ".dat");
	strcat(cpixel_file3,     ".xyz");
	strcat(cpixel_file4,     ".dat");
	strcat(cpixel_file5,     ".dat");
	strcat(cpixel_file6,     ".xyz");

	ofstream out_slab_in(cpixel_file2);
	ofstream inside_neu(cpixel_file3);
	ofstream out_slab_notinterst(cpixel_file4);
	ofstream out_slab_interst(cpixel_file5);
	ofstream out_slab_bombloc(cpixel_file6);

	float x_trans;
	for (i=0; i<=neural_array.last_one; i++)
	{
	  // WRITE NEURON LOCATIONS THAT ARE INSIDE THE SLAB
	  if (neural_array[i].isAlive() && neural_array[i].isInside())
	  {
		// WRITE 3D POSITIONS
		if (dq.write_files)
		//	inside_neu << neural_array[i].getWrapPosition(dq.L_box) << "\n";
			inside_neu << neural_array[i].getRotatedPosition() << "\n";

		// PROJECT TO THE X'-Z' PLANE
		x_trans = (neural_array[i].getRotatedPosition().x - dq.L_box/2.0) * cos(dq.angle) 
		  + (neural_array[i].getRotatedPosition().y - dq.L_box/2.0) * sin(dq.angle)
		  + dq.L_tissue/2.0;

		// WRITE 2D POSITIONS FOR ALL NEURONS
		out_slab_in << x_trans << " " << neural_array[i].getRotatedPosition().z << " " << neural_array[i].getMarked() << "\n";
		// WRITE ONLY INTERSTITIAL OR NOT INTERSTITIAL
		if (dq.write_files && !neural_array[i].isInterst())
		{
		  	// WRITE COLUMNAR AND NOT INTERSTITIAL
			out_slab_notinterst << x_trans << " " << neural_array[i].getRotatedPosition().z << "\n";
		}
		if (dq.write_files &&  neural_array[i].isInterst())
		{
		  	// WRITE INTERSTITIAL
			out_slab_interst << x_trans << " " << neural_array[i].getRotatedPosition().z << "\n";
		}
	  }

	  // WRITE BOMB LOCATIONS THAT ARE INSIDE THE SLAB
	  // 		NOTE THAT FOR BOMB LOCATIONS,  neural_array[i].isAlive()==0
	  if (neural_array[i].isBombLoc() && neural_array[i].isInside())
	  {
		// WRITE 3D POSITIONS
		if (dq.write_files)
		//	out_slab_bombloc << neural_array[i].getWrapPosition(dq.L_box) << "\n";
			out_slab_bombloc << neural_array[i].getRotatedPosition() << "\n";
	  }
	}
	inside_neu.close();
	out_slab_in.close();
	out_slab_notinterst.close();
	out_slab_interst.close();
	out_slab_bombloc.close();
	delete [] cpixel_file2;
	delete [] cpixel_file3;
	delete [] cpixel_file4;
	delete [] cpixel_file5;
	delete [] cpixel_file6;
	//}
	//
	// B. NOW, WRITE THE 3D BLOCK GLIA INSIDE THE SLAB TO FILE.
	// ALSO, TRANSFORM FROM THE X-Y-Z BASIS TO THE X'-Y'-Z' ROTATED BASIS.
	// GENERATE X'-Z' DATA ONLY (TISSUE COLLAPSE) WHERE Z=Z'.
	if (dq.put_glia == 1)
	{
	  cout << " [Projecting glia and writing out.]\n";
	  //if (dq.write_files)
	  //{
	  char pre_file7[40]     = "slab_2DGlia_";
	  char pre_file8[40]     = "slab_3DGlia_";
	  strcat(pre_file7, file_lab);
	  strcat(pre_file7, "_");
	  strcat(pre_file8, file_lab);
	  strcat(pre_file8, "_");
	  char *cpixel_file7     = NULL;
	  cpixel_file7     = conc_char_2_int(pre_file7,dq.num_steps,3);
	  char *cpixel_file8     = NULL;
	  cpixel_file8     = conc_char_2_int(pre_file8,dq.num_steps,3);
	  strcat(cpixel_file7,     ".dat");
	  strcat(cpixel_file8,     ".xyz");

	  ofstream out_slab_in_glia(cpixel_file7);
	  ofstream inside_glia(cpixel_file8);

	  float x_trans;
	  for (i=0; i<=glial_array.last_one; i++)
	  {
		// WRITE GLIA LOCATIONS THAT ARE INSIDE THE SLAB
		if (glial_array[i].isAlive() && glial_array[i].isInside())
		{
		  // WRITE 3D POSITIONS
		  if (dq.write_files)
			//	inside_glia << glial_array[i].getWrapPosition(dq.L_box) << "\n";
			inside_glia << glial_array[i].getRotatedPosition() << "\n";

		  // PROJECT TO THE X'-Z' PLANE
		  x_trans = (glial_array[i].getRotatedPosition().x - dq.L_box/2.0) * cos(dq.angle) 
			+ (glial_array[i].getRotatedPosition().y - dq.L_box/2.0) * sin(dq.angle)
			+ dq.L_tissue/2.0;

		  // WRITE 2D POSITIONS FOR ALL GLIA
		  out_slab_in_glia << x_trans << " " << glial_array[i].getRotatedPosition().z << "\n";
		}
	  }
	  inside_glia.close();
	  out_slab_in_glia.close();
	  delete [] cpixel_file7;
	  delete [] cpixel_file8;
	  //}
	}
	// END: WRITE THE 3D BLOCK NEURONS AND GLIA INSIDE THE SLAB TO FILE.


	// ==========================================================
	// WRITE THE 2D PROJECTION ALSO IN PIC FORMAT 
	if (dq.write_files)
	{
	  // 1. WRITE NEURON DATA
	  char pic_file[40]     = "slab_2D_";
	  strcat(pic_file, file_lab);
	  strcat(pic_file, "_");
	  char *cpixel_pic     = NULL;
	  cpixel_pic     = conc_char_2_int(pic_file,dq.num_steps,3);
	  strcat(cpixel_pic,     ".pic");

	  // IMAGE SIZE FOR THE PIC FILE
	  int image_size  = 512;
	  // SIZE OF THE PARTICLE IN THE PIC FILE
	  //float part_size = 3.;
	  //float part_size = 5.;
	  float part_size = dq.neu_radius;

	  // WRITE PIC FILE AS ONE SLICE
	  writePIC_neurons_cb(neural_array, dq, cpixel_pic, image_size, part_size);
	  // WRITE PIC FILE AS TWO SLICES, THE SECOND ONE BLANK
	  //write_twoPICs_cb(neural_array, dq, cpixel_pic, image_size, part_size);
	  delete [] cpixel_pic;

	  // 2. NOW, WRITE GLIA DATA
	  if (dq.put_glia == 1)
	  {
		char pic_file_glia[40]     = "slab_2DGlia_";
		strcat(pic_file_glia, file_lab);
		strcat(pic_file_glia, "_");
		char *cpixel_pic_glia     = NULL;
		cpixel_pic_glia     = conc_char_2_int(pic_file_glia,dq.num_steps,3);
		strcat(cpixel_pic_glia,     ".pic");

		// IMAGE SIZE FOR THE PIC FILE
		int image_size_glia  = 512;
		// SIZE OF THE PARTICLE IN THE PIC FILE
		float part_size_glia = dq.glia_radius;

		// WRITE PIC FILE AS ONE SLICE
		writePIC_glia_cb(glial_array, dq, cpixel_pic_glia, image_size_glia, part_size_glia);
		// WRITE PIC FILE AS TWO SLICES, THE SECOND ONE BLANK
		// -- DOES NOT EXIST YET
		//write_twoPICs_glia_cb(glial_array, dq, cpixel_pic_glia, image_size_glia, part_size_glia);
	  	delete [] cpixel_pic_glia;
	  }
	}

	// ===============================================
	// DEALLOCATE MEMORY
	delete [] x_cp;
	delete [] y_cp;

	for (i=0; i<num_eveneu_x; i++)
	{
	  delete [] basemat_even_x[i];
	  delete [] basemat_even_y[i];
	}
	delete [] basemat_even_x;
	delete [] basemat_even_y;

	for (i=0; i<num_oddneu_x; i++)
	{
	  delete [] basemat_odd_x[i];
	  delete [] basemat_odd_y[i];
	}
	delete [] basemat_odd_x;
	delete [] basemat_odd_y;


	delete [] file_lab;
	delete [] the_time;

	param_out.close();
	delete [] param_file;

  } // END OF CONFIGURATION LOOP WITH COUNTER cf

#ifdef __MPI__
    MPI_Finalize();
#endif

  return 0;
  // RETURN RANDOM LONG INTEGER. USED WHEN GENERATING A SEQUENCE
  // OF CONFIGURATIONS.
  //return -rand_long();
}

void mean_neural_displacement(
	double& ave_Rdis, double& stdev, 
	double& ave_Rdis_XY, double& stdev_XY, 
	double& ave_Rdis_Z, double& stdev_Z, 
	infobjarray<neuron,NEU_ARRAY>& neural_array, 
	data& dq, float box)
  // RETURN AVERAGE AND STANDARD DEVIATION.
  // CALCULATE DISPLACEMENT BY AVERAGING THE SQUARE OF THE DISPLACEMENT
  // AND THEN CALCULATING THE SQUARE ROOT.
  // THIS IS THE "SPREAD" IN A CLASSICAL RANDOM WALK.
{
  ave_Rdis    = 0.;
  ave_Rdis_XY = 0.;
  ave_Rdis_Z  = 0.;
  stdev = 0.;
  double tmp_disp2, tmp_disp_XY2, tmp_disp_Z2;
  float boxinv = 1.0 / box;
  int i,j;

  for (i=0; i<=neural_array.last_one; i++)
  {
	if (neural_array[i].isAlive() )
	{
	  // 3D DISPLACEMENT
	  //ave_Rdis += neural_array[i].getImageDisplacement(boxinv, box);
	  tmp_disp2 = neural_array[i].getDisplacement2();
	  ave_Rdis += tmp_disp2;
	  stdev    += CUAD( tmp_disp2 );

	  // 2D - XY DISPLACEMENT
	  tmp_disp_XY2  = neural_array[i].getXYDisplacement2();
	  ave_Rdis_XY += tmp_disp_XY2;
	  stdev_XY    += CUAD( tmp_disp_XY2 );

	  // 1D - Z DISPLACEMENT
	  tmp_disp_Z2 = neural_array[i].getZDisplacement2();
	  ave_Rdis_Z += tmp_disp_Z2;
	  stdev_Z    += CUAD( tmp_disp_Z2 );
	}
  }

  // 3D DISPLACEMENT
  double dummy = 0.000001;
  if ( ((stdev-ave_Rdis)>0.) && (ABS(stdev-ave_Rdis)>dummy) )
  {
	ave_Rdis /= double(dq.num_neu);
	stdev    /= double(dq.num_neu);
	stdev     = sqrt( stdev - CUAD(ave_Rdis) );
	// N-1 STDDEV
	stdev     = stdev * sqrt(  dq.num_neu /(dq.num_neu-1) );

	// CALCULATE SPREAD
	ave_Rdis = sqrt( ave_Rdis );
	//	  not sure of SQRT of stdev is correct.
	stdev    = sqrt( stdev    );
  }
  else
  {
	ave_Rdis /= double(dq.num_neu);
  	stdev     = 0.;

	// CALCULATE SPREAD
	ave_Rdis = sqrt( ave_Rdis );
  }

  // 2D - XY DISPLACEMENT
  dummy = 0.000001;
  if ( ((stdev_XY-ave_Rdis_XY)>0.) && (ABS(stdev_XY-ave_Rdis_XY)>dummy) )
  {
	ave_Rdis_XY /= double(dq.num_neu);
	stdev_XY    /= double(dq.num_neu);
	stdev_XY     = sqrt( stdev_XY - CUAD(ave_Rdis_XY) );
	// N-1 STDDEV
	stdev_XY     = stdev_XY * sqrt( dq.num_neu /(dq.num_neu-1) );

	// CALCULATE SPREAD
	ave_Rdis_XY = sqrt( ave_Rdis_XY );
	//	  not sure of SQRT of stdev is correct.
	stdev_XY    = sqrt( stdev_XY    );
  }
  else
  {
	ave_Rdis_XY /= double(dq.num_neu);
  	stdev_XY     = 0.;

	// CALCULATE SPREAD
	ave_Rdis_XY = sqrt( ave_Rdis_XY );
  }

  // 1D - Z DISPLACEMENT
  dummy = 0.000001;
  if ( ((stdev_Z-ave_Rdis_Z)>0.) && (ABS(stdev_Z-ave_Rdis_Z)>dummy) )
  {
	ave_Rdis_Z /= double(dq.num_neu);
	stdev_Z    /= double(dq.num_neu);
	stdev_Z     = sqrt( stdev_Z - CUAD(ave_Rdis_Z) );
	// N-1 STDDEV
	stdev_Z     = stdev_Z * sqrt( dq.num_neu /(dq.num_neu-1) );

	// CALCULATE SPREAD
	ave_Rdis_Z = sqrt( ave_Rdis_Z );
	//	  not sure of SQRT of stdev is correct.
	stdev_Z    = sqrt( stdev_Z    );
  }
  else
  {
	ave_Rdis_Z /= double(dq.num_neu);
  	stdev_Z     = 0.;

	// CALCULATE SPREAD
	ave_Rdis_Z = sqrt( ave_Rdis_Z );
  }
}

void write_displacement_distribution(infobjarray<neuron,NEU_ARRAY>& neural_array, char* the_file)
  // WRITES TO FILE THE DISTRIBUTION OF DISPLACEMENTS
{
  int i,j;

  ofstream out_file(the_file);

  for (i=0; i<=neural_array.last_one; i++)
  {
	if (neural_array[i].isAlive() )
		out_file << neural_array[i].getDisplacement() << "\n";
  }
  out_file.close();
}

void gen_rand_angles(data& dq)
{
  // VERTICAL ANGLE OF THE SLAB, IN THE XY PLANE.
  // THE ANGLE IS IN THE XZ PLANE, AND IS EQUIVALENT
  // TO THE ANGLE OF THE BLADE CUTTING THE SLICE.
  if (dq.angle<0.0)
  {
	dq.angle = random_angle(0.,180.0);
	//cout << "Running with random angle: " << dq.angle << "\n";
  }
  dq.angle = RAD(dq.angle);

  // AZIMUTH ANGLE FOR SIDEWAYS CUT
  if (dq.angle_phi<0.0)
  {
	dq.angle_phi = random_angle(dq.angle_rndphi_low, dq.angle_rndphi_high);
  }
  dq.angle_phi = RAD(dq.angle_phi);

  // IN-PLANE ANGLE FOR SIDEWAYS CUT
  if (dq.angle_theta<0.0)
  {
	dq.angle_theta = random_angle(0., 360.0);
	//cout << "Running with random theta angle: " << dq.angle_theta << "\n";
  }
  dq.angle_theta = RAD(dq.angle_theta);
}

void write_block_to_file(infobjarray<neuron,NEU_ARRAY>& neural_array, data& dq, char* file_n)
{
  int i;
  ofstream init_neu(file_n);
  for (i=0; i<=neural_array.last_one; i++)
  {
	// WRITE TO FILE THE WRAPPED POSITION INSIDE THE BOX.
	if (neural_array[i].isAlive())
		init_neu << neural_array[i].getWrapPosition(dq.L_box) << "\n";
  }
  init_neu.close();
}

void write_block_to_file_colID(infobjarray<neuron,NEU_ARRAY>& neural_array, data& dq, char* file_n)
{
  int i;
  ofstream init_neu(file_n);
  for (i=0; i<=neural_array.last_one; i++)
  {
	// WRITE TO FILE THE WRAPPED POSITION INSIDE THE BOX.
	if (neural_array[i].isAlive())
		init_neu 
		  << neural_array[i].getMarked() << " "
		  << neural_array[i].getWrapPosition(dq.L_box) << "\n";
  }
  init_neu.close();
}

int neuralForcesZero(infobjarray<neuron,NEU_ARRAY>& neural_array)
{
  int i,j;
  int forces_zero = 1;

  vec zero(0.,0.,0.);

  for (i=0; i<=neural_array.last_one; i++)
  {
	if ( neural_array[i].isAlive() && ((neural_array[i].forces) != zero) )
	  forces_zero = 0;
  }

  return forces_zero;
}

void updateOriginalPositions(infobjarray<neuron,NEU_ARRAY>& neural_array)
{
  for (int i=0; i<=neural_array.last_one; i++)
	if (neural_array[i].isAlive())
		neural_array[i].updateOriginalPosition();
}

void write_density_distribution(infobjarray<neuron,NEU_ARRAY>& neural_array, data& dq, char* file_n)
  // DIVIDES THE BOX INTO CELLS AND WRITES THE CELL NEURONAL DENSITIES
{
  int i,j,k;

  double cell_size = dq.dens_cell_size;
  
  double vol_cell  = cell_size * cell_size * cell_size;

  // NUMBER OF CELLS
  int Nc = int(dq.L_box/cell_size) + 1;

  // CREATE MATRIX
  int ***densmap = new int**[Nc];
  for (i=0; i<Nc; i++)
  {
	densmap[i] = new int*[Nc];
	for (j=0; j<Nc; j++)
	  densmap[i][j] = new int[Nc];
  }
  // INITIALIZE
  for (i=0; i<Nc; i++)
	for (j=0; j<Nc; j++)
	  for (k=0; k<Nc; k++)
		densmap[i][j][k] = 0;

  // LOOP OVER EACH NEURON AND ADD TO CELL
  for (i=0; i<=neural_array.last_one; i++)
  {
	// WRITE TO FILE THE WRAPPED POSITION INSIDE THE BOX.
	if (neural_array[i].isAlive())
	{
	  vec pos = neural_array[i].getWrapPosition(dq.L_box);
	  int ni = int(pos.x/cell_size);
	  int nj = int(pos.y/cell_size);
	  int nk = int(pos.z/cell_size);

	  // ADD TO DENSITY MATRIX
	  densmap[ni][nj][nk]++;
	}
  }

  // WRITE TO FILE
  ofstream densmap_out(file_n);
  densmap_out.precision(10);
  //   DO NOT WRITE THE EDGE CELLS, THAT HAVE ARTIFICIALLY LOW DENSITY
  //   BECAUSE L_box MAY NOT BE A MULTIPLE OF cell_size.
  for (i=0; i<Nc-1; i++)
	for (j=0; j<Nc-1; j++)
	  for (k=0; k<Nc-1; k++)
	  {
		densmap_out << i << " "
		            << j << " "
		            << k << " "
					<< densmap[i][j][k] / vol_cell << "\n";
	  }
  densmap_out.close();


  // DEALLOCATE
  for (i=0; i<Nc; i++)
  {
	for (j=0; j<Nc; j++)
	{
	  delete [] densmap[i][j];
	}
	delete [] densmap[i];
  }
  delete [] densmap;

}

/*
 * -------------------------------------------------------
 */
