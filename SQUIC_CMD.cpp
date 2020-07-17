// Data structures, algos and numerics
#include <string>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <ctime>
#include <math.h>
#include <iomanip>

#ifdef _WIN32
#include <Windows.h>
#else
#include <unistd.h>
#endif


#define IGNORE_STRING	"NULL"

#define FILE_RETRY_SEC	2


extern "C" {

	void SQUIC_C(
	    int p, int n, double* Y_ptr,
	    double lambda,
	    int* LambdaMatrix_ij, double* LambdaMatrix_val, int LambdaMatrix_nnz,
	    int max_iter, double drop_tol, double term_tol, int verbose,
	    int* X0_ij,            double* X0_val,            int X0_nnz,
	    int* W0_ij,            double* W0_val,            int W0_nnz,
	    int*&  X_ij,            double*&  X_val,            int& X_nnz,
	    int*&  W_ij,            double*&  W_val,            int& W_nnz,
	    int&    stat_niters,
	    double* stat_times, // length of 7: [total,cov,itr,chol,inv,lns,upd]
	    double& stat_dgap,
	    double*& stat_obj
	);
}

/**
 * @brief      Return the number of nnz in COO matrix
 *
 * @param[in]  filename  The filename of COO format matrix.
 *
 * @return     { int: Number of nonzeros }
 */
int read_COO_nnz(std::string filename) {

	int count = 0;
	std::ifstream inputfile(filename);


	if (!inputfile.is_open()) {
		std::cerr << "[Warning] Unable to open file " << filename << std::endl;
		std::cerr << "[Warning] Sleeping for " << FILE_RETRY_SEC << "Seconds" << std::endl;
		sleep(FILE_RETRY_SEC);
		if (!inputfile.is_open()) {
			std::cerr << "[Error] Unable to open file " << filename << std::endl;
			exit(1);   // call system to stop
		}
	}

	// Count the number of lines = nnz
	std::string line;
	while (std::getline(inputfile, line)) {
		++count;
	}

	inputfile.close();

	return count;
};


/**
 * @brief      Reads a coo matrix.
 *
 * @param[in]  filename     The filename
 * @param      ij           Index array
 * @param      val          Value array
 * @param[in]  base_offset  Index basis offset
 */
void read_COO_Mat(std::string filename, int* ij, double* val, int base_offset) {

	int i, j, k;
	double v;

	std::ifstream inputfile(filename);
	if (!inputfile.is_open()) {
		std::cerr << "[Warning] Unable to open file " << filename << std::endl;
		std::cerr << "[Warning] Sleeping for " << FILE_RETRY_SEC << "Seconds" << std::endl;
		sleep(FILE_RETRY_SEC);
		if (!inputfile.is_open()) {
			std::cerr << "[Error] Unable to open file " << filename << std::endl;
			exit(1);   // call system to stop
		}
	}

	k = 0;
	while (inputfile >> i >> j >> v) {
		ij[2 * k]     = i - base_offset;
		ij[2 * k + 1] = j - base_offset;
		val[k]        = v;
		++k;
	};

	inputfile.close();
};


/**
 * @brief      Writes a coo matrix.
 *
 * @param[in]  filename     Filename of COO Matrix
 * @param      X_ij         Index array
 * @param      X_val        Value array
 * @param[in]  nnz          Number of nonzeros
 * @param[in]  base_offset  Index basis offset
 */
void write_COO_Mat(std::string filename, int* ij, double* val, int nnz, int base_offset) {

	// We return it the same way it was read.
	// The ouput of SQUIC is zero based so switch to what the input was

	std::ofstream outputfile(filename);
	if (!outputfile.is_open()) {
		std::cerr << "[Warning] Unable to open file " << filename << std::endl;
		std::cerr << "[Warning] Sleeping for " << FILE_RETRY_SEC << "Seconds" << std::endl;
		sleep(FILE_RETRY_SEC);
		if (!outputfile.is_open()) {
			std::cerr << "[Error] Unable to open file " << filename << std::endl;
			exit(1);   // call system to stop
		}
	}


	for (int i = 0; i < nnz; ++i) {
		outputfile << ij[2 * i] + base_offset  << "	" << ij[2 * i + 1] + base_offset  << "	" << val[i] << std::endl;
	};
	outputfile.close();
};



int main(int argc, char* argv[]) {


	if (argc != 15 + 1) {
		std::cerr << "SQUIC_CMD : p n Y_loc lambda LambdaMatrix_loc max_iter drop_tol term_tol X0_loc W0_loc index_offset verbose X_loc W_loc " << std::endl;
		std::cerr << "[integer:p] Number of random variables" << std::endl;
		std::cerr << "[integer:n] Number of samples." << std::endl;
		std::cerr << "[string:Y_loc] Input loc/of/file/Y.dat" << std::endl;
		std::cerr << "[double:lambda] Scalar lambda value." << std::endl;
		std::cerr << "[string:LambdaMatrix_loc] loc/of/file/LambdaMatrix.dat" << std::endl;
		std::cerr << "[integer:max_iter] Max iteration." << std::endl;
		std::cerr << "[double:drop_tol] Drop out tolderance." << std::endl;
		std::cerr << "[double:term_tol] Terminal Tolerence" << std::endl;
		std::cerr << "[string:X0_loc] Input loc/of/file/X0.dat" << std::endl;
		std::cerr << "[string:W0_loc] Input loc/of/file/W0.dat" << std::endl;
		std::cerr << "[integer:index_offset] Offset of indexing" << std::endl;
		std::cerr << "[integer:verbose] Verbosity" << std::endl;
		std::cerr << "[string:X_loc] Output loc/of/file/X.dat" << std::endl;
		std::cerr << "[string:W_loc] Output loc/of/file/W.dat" << std::endl;
		std::cerr << "[string:log_loc] Output loc/of/file/log.dat" << std::endl;
		std::cerr << "NOTE: ~:Optional input. Can be ignored by using (Case Sensiteve): NULL" << std::endl;
		exit(1);
	}

	int p							= atoi(argv[1]);
	int n							= atoi(argv[2]);
	std::string Y_loc				= argv[3];
	double lambda					= atof(argv[4]);
	std::string LambdaMatrix_loc	= argv[5];
	int max_iter					= atoi(argv[6]);
	double drop_tol					= atof(argv[7]);
	double term_tol					= atof(argv[8]);
	std::string X0_loc				= argv[9];
	std::string W0_loc				= argv[10];
	int index_offset				= atoi(argv[11]);
	int verbose						= atoi(argv[12]);
	std::string X_loc				= argv[13];
	std::string W_loc				= argv[14];
	std::string log_loc				= argv[15];

	// Read Data matrix
	std::fstream Y_file(Y_loc, std::ios_base::in);
	if (!Y_file.is_open()) {
		std::cerr << "[Warning] Unable to open file " << Y_loc << std::endl;
		std::cerr << "[Warning] Sleeping for " << FILE_RETRY_SEC << "Seconds" << std::endl;
		sleep(FILE_RETRY_SEC);
		if (!Y_file.is_open()) {
			std::cerr << "[Error] Unable to open file " << Y_loc << std::endl;
			exit(1);   // call system to stop
		}
	}
	double* Y = new double[p * n];
	int i = 0;
	double Y_val;
	while (Y_file >> Y_val) {
		Y[i] = Y_val;
		i++;
	}
	Y_file.close();


	int* X0_ij;
	double* X0_val;
	int X0_nnz;

	int* W0_ij;
	double* W0_val;
	int W0_nnz;

	int* LambdaMatrix_ij;
	double* LambdaMatrix_val;
	int LambdaMatrix_nnz;

	int* X_ij;
	double* X_val;
	int X_nnz;

	int* W_ij;
	double* W_val;
	int W_nnz;

	int* IGNORE_I;
	double* IGNORE_D;


	// stat variables
	int     stat_niters;
	double  stat_times[7];
	double  stat_dgap;
	double* stat_obj;        // size of stat_niters

	bool use_scalar_lambda = LambdaMatrix_loc == IGNORE_STRING;
	bool use_identiy_intial_guess = X0_loc == IGNORE_STRING || W0_loc == IGNORE_STRING;


	// Using scalar lambda & identity intial guess
	if (use_scalar_lambda && use_identiy_intial_guess ) {
		SQUIC_C(
		    p, n, Y,
		    lambda,
		    IGNORE_I,	IGNORE_D,	-1,
		    max_iter,	drop_tol,	term_tol,	verbose,
		    IGNORE_I,	IGNORE_D,	-1,
		    IGNORE_I,	IGNORE_D,	-1,
		    X_ij,		X_val,		X_nnz,
		    W_ij,		W_val,		W_nnz,
		    stat_niters,
		    stat_times, // length of 7: [total,cov,itr,chol,inv,lns,upd]
		    stat_dgap,
		    stat_obj
		);
	}
	// Using Matrix lambda & identity intial guess
	else if (!use_scalar_lambda && use_identiy_intial_guess) {

		LambdaMatrix_nnz  = read_COO_nnz(LambdaMatrix_loc);
		LambdaMatrix_ij   = (int*)    malloc(LambdaMatrix_nnz * 2 * sizeof(int));
		LambdaMatrix_val  = (double*)  malloc(LambdaMatrix_nnz *    sizeof(double));
		read_COO_Mat(LambdaMatrix_loc, LambdaMatrix_ij, LambdaMatrix_val, index_offset);

		SQUIC_C(
		    p, n, Y,
		    lambda,
		    LambdaMatrix_ij,	LambdaMatrix_val, 	LambdaMatrix_nnz,
		    max_iter, drop_tol, term_tol, verbose,
		    IGNORE_I,           IGNORE_D,			-1,
		    IGNORE_I,           IGNORE_D,			-1,
		    X_ij,				X_val,				X_nnz,
		    W_ij,				W_val,				W_nnz,
		    stat_niters,
		    stat_times, // length of 7: [total,cov,itr,chol,inv,lns,upd]
		    stat_dgap,
		    stat_obj
		);

		delete[] LambdaMatrix_ij;
		delete[] LambdaMatrix_val;
	}

	// Using Matrix lambda & provided intial guess
	else if (!use_scalar_lambda && !use_identiy_intial_guess) {

		LambdaMatrix_nnz  = read_COO_nnz(LambdaMatrix_loc);
		LambdaMatrix_ij   = (int*)    malloc(LambdaMatrix_nnz * 2 * sizeof(int));
		LambdaMatrix_val  = (double*)  malloc(LambdaMatrix_nnz *    sizeof(double));
		read_COO_Mat(LambdaMatrix_loc, LambdaMatrix_ij, LambdaMatrix_val, index_offset);

		X0_nnz  = read_COO_nnz(X0_loc);
		X0_ij   = (int*)    malloc(X0_nnz * 2 * sizeof(int));
		X0_val   = (double*)  malloc(X0_nnz *    sizeof(double));
		read_COO_Mat(X0_loc, X0_ij, X0_val, index_offset);

		W0_nnz  = read_COO_nnz(W0_loc);
		W0_ij   = (int*)    malloc(W0_nnz * 2 * sizeof(int));
		W0_val  = (double*)  malloc(W0_nnz *    sizeof(double));
		read_COO_Mat(W0_loc, W0_ij, W0_val, index_offset);

		SQUIC_C(
		    p, n, Y,
		    lambda,
		    LambdaMatrix_ij,	LambdaMatrix_val, 	LambdaMatrix_nnz,
		    max_iter, drop_tol, term_tol, verbose,
		    X0_ij,          	X0_val,				X0_nnz,
		    W0_ij,           	W0_val,				W0_nnz,
		    X_ij,				X_val,				X_nnz,
		    W_ij,				W_val,				W_nnz,
		    stat_niters,
		    stat_times, // length of 7: [total,cov,itr,chol,inv,lns,upd]
		    stat_dgap,
		    stat_obj
		);

		delete[] LambdaMatrix_ij;
		delete[] LambdaMatrix_val;

		delete[] X0_ij;
		delete[] X0_val;

		delete[] W0_ij;
		delete[] W0_val;
	}

	// Using scalar lambda & provided intial guess
	else if (use_scalar_lambda && !use_identiy_intial_guess) {

		X0_nnz  = read_COO_nnz(X0_loc);
		X0_ij   = (int*)    malloc(X0_nnz * 2 * sizeof(int));
		X0_val   = (double*)  malloc(X0_nnz *    sizeof(double));
		read_COO_Mat(X0_loc, X0_ij, X0_val, index_offset);

		W0_nnz  = read_COO_nnz(W0_loc);
		W0_ij   = (int*)    malloc(W0_nnz * 2 * sizeof(int));
		W0_val  = (double*)  malloc(W0_nnz *    sizeof(double));
		read_COO_Mat(W0_loc, W0_ij, W0_val, index_offset);


		SQUIC_C(
		    p, n, Y,
		    lambda,
		    IGNORE_I,	IGNORE_D, 	-1,
		    max_iter, drop_tol, term_tol, verbose,
		    X0_ij,          	X0_val,				X0_nnz,
		    W0_ij,           	W0_val,				W0_nnz,
		    X_ij,				X_val,				X_nnz,
		    W_ij,				W_val,				W_nnz,
		    stat_niters,
		    stat_times, // length of 7: [total,cov,itr,chol,inv,lns,upd]
		    stat_dgap,
		    stat_obj
		);

		delete[] X0_ij;
		delete[] X0_val;

		delete[] W0_ij;
		delete[] W0_val;

	} else {
		std::cerr << "Something went wrong with the options..." << std::endl;
		exit(0);
	}


//	std::cout << "X_nnz_sym:" << X_nnz << std::endl;
//	std::cout  << "W_nnz_sym:" << W_nnz << std::endl;


//	for (int i = 0; i < X_nnz; ++i) {
//		std::cout << X_ij[2 * i] << " " << X_ij[2 * i + 1] << " " << X_val[i] << std::endl;
//	}

	write_COO_Mat(X_loc, X_ij, X_val, X_nnz, index_offset);
	write_COO_Mat(W_loc, W_ij, W_val, W_nnz, index_offset);

	std::ofstream log_file(log_loc,    std::ios::out | std::ios::trunc);
	log_file << "{" << std::endl;
	log_file << "\"time_tamp\":" << std::time(0) << "," << std::endl;
	log_file << "\"p\":" << p << "," << std::endl;
	log_file << "\"n\":" << n << "," << std::endl;
	log_file << "\"lambda\":" << lambda << "," << std::endl;
	log_file << "\"LambdaMatrix_loc\":\"" << LambdaMatrix_loc << "\"" << "," << std::endl;
	log_file << "\"max_iter\":" << max_iter << "," << std::endl;
	log_file << "\"drop_tol\":" << drop_tol << "," << std::endl;
	log_file << "\"term_tol\":" << term_tol << "," << std::endl;
	log_file << "\"verbose\":" << verbose << "," << std::endl;
	log_file << "\"X0_loc\":\"" << X0_loc << "\"" << "," << std::endl;
	log_file << "\"W0_loc\":\"" << W0_loc << "\"" << "," << std::endl;
	log_file << "\"X_nnz\":" << X_nnz << "," << std::endl;
	log_file << "\"W_nnz\":" << W_nnz << "," << std::endl;
	log_file << "\"niters\":" << stat_niters << "," << std::endl;
	log_file << "\"time_total\":" << stat_times[0] << "," << std::endl;
	log_file << "\"time_cov\":" << stat_times[1] << "," << std::endl;
	log_file << "\"time_itr\":" << stat_times[2] << "," << std::endl;
	log_file << "\"time_chol\":" << stat_times[3] << "," << std::endl;
	log_file << "\"time_inv\":" << stat_times[4] << "," << std::endl;
	log_file << "\"time_lns\":" << stat_times[5] << "," << std::endl;
	log_file << "\"time_upd\":" << stat_times[6] << "," << std::endl;
	log_file << "\"dgap\":" << stat_dgap << "," << std::endl;
	log_file << "\"obj\":[";
	for (int i = 0; i < stat_niters - 1; ++i) {
		log_file << std::setprecision(12) << stat_obj[i] << ",";
	}
	log_file << stat_obj[stat_niters - 1] << "]" << "," << std::endl;

	double X_norm_1 = 0;
	for (int i = 0; i < X_nnz; ++i) {
		X_norm_1 += fabs(X_val[i]);
	}
	log_file << "\"X_norm_1\":" << X_norm_1 << "," << std::endl;

	double W_norm_1 = 0;
	for (int i = 0; i < W_nnz; ++i) {
		W_norm_1 += fabs(W_val[i]);
	}
	log_file << "\"W_norm_1\":" << W_norm_1 << std::endl;
	log_file << "}" << std::endl;

	delete[] X_ij;
	delete[] X_val;

	delete[] W_ij;
	delete[] W_val;

	delete[] Y;

	return 0;
}
