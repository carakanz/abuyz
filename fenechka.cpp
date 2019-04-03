// June13_BlockPxPy.cpp : Defines the entry point for the console application.
//
//


#include "cnpy.h"
#include <cstdlib>
#include "iostream" //»менно iostream, а не iostream.h
#include "fstream"
#include "math.h"
#include "complex"
#include <cmath>
#include <iomanip>
#include <time.h>
#include <map>
#include <vector>
#include <string.h>
#include <limits.h>


//#include <thread>
//#pragma STDC CX_LIMITED_RANGE on/off

#define PI 3.1415926536f
#define sqrt2 1.41421356f

#define hbar 1.055e-27f
#define ee 4.803e-10f
#define cc 2.998e10f
#define eps0 8.5f
#define mm 2.733e-28f  //me=0.3 m0  in ZnO
#define muBohr 9.274e-21f
#define gfact 2.f
#define meV_erg 1.602e-15f

#define float_err 1e-5f



using namespace std;

typedef unsigned int uint;
typedef unsigned long long ullong;
//////////////////////////////////////////////////////////
//////F U N C T I O N S   D E C L A R A T I O N S /////////////////////////////
//////////////////////////////////////////////////////////
//////F U N C T I O N S   D E C L A R A T I O N S /////////////////////////////


//ullong* Construct_Fermi_Basis_L(uint Ne, uint NS, uint NLL, float Sz);
void Construct_Fermi_Basis_L(std::vector<ullong>& s_Sz, uint Ne, uint NS, uint NLL, float Sz);
void display_vect(ullong vv, uint NS, uint NLL);
ullong Birth(ullong vect, uint pos);
int SignBirth(ullong vector, uint pos);
ullong Death(ullong vect, uint pos);
int SignDeath(ullong vector, uint pos);
ullong TranslVectorSigned(ullong vector, int &ssign, uint NS, uint NLL);//here we perform a collinear translation with the sign calculation
uint Calcul_qy_for_vector(ullong vector, uint Ne, uint NS, uint NLL);

uint Filter_Basis_qy(const std::vector<ullong>& in_basSz, std::vector<uint>& out_blockBasis, std::vector<int>& out_signs, std::vector<uint>& out_sizes,
	uint Ne, uint NU, uint NS, uint NLL, float Sz, uint qy); //create sub_basis of states with the fixed qy-projection
//uint Filter_Basis_qy(ullong* in_basSz, uint* out_blockBasis, int* out_signs, uint* out_sizes,
//	uint Ne, uint NU, uint NS, uint NLL, float Sz, uint qy); //create sub_basis of states with the fixed qy-projection

uint Basis_size_any_NU_Sz(uint _Ne, uint _NS, uint _NLL, float Sz);
uint Combinat_N_m(uint N, uint m);
uint index_vector(ullong vect, const std::vector<ullong>& arr_vect, uint NMAX);

uint DeltaSymbol(int pos1, int pos2, uint Nphi);
float FFnmReal(float x, uint n, uint m);
ullong TwoParticleMixingOperator(ullong vector, int &ssign, uint i1, uint i2, uint i3, uint i4, uint n1, uint n2, uint n3, uint n4, uint spin14, uint spin23, uint Nphi);
complex<float> MatrixElement(int i1, int i2, int i3, int i4,
	uint n1, uint n2, uint n3, uint n4,
	uint Nphi, float Lb, float LX, float LY, float ECl, const std::vector<float>& eps_s, const std::vector<float>& FormFactArray, uint Q_cutoff);

//float* calc_eps_s_full(uint Nphi, float Lb, float LX, float LY, float ECl_ECr_ratio, uint _NLL);
//APPLICABLE ONLY FOR EVEN NUMBER OF ELECTRONS!!!  :
float Sz_full_represent(ullong vector, uint _Nphi, uint _NLL); // it means that both spin states of each LL are considered
float CalcSingleParticle_En(ullong vector, uint _NLL, uint _NS, float Ecr, float Ezm);
//float eps_scr_symmetric(float q, float ratioE, uint _NLL); //just analitic function

uint factorial(uint N);
float sq(float x);
uint max(uint a, uint b);
uint min(uint a, uint b);
int sign_u(uint k, uint m);  //returns sign of (k-m)
complex<float> pow_cmplxint(complex<float> z, uint t);
complex<float> conjugate(complex<float> z); //returns complex conjugate
void round_cmplx(complex<float> &z);

//SCREENING functions

float FormFact(float q, float w_lb);
void calc_FormFactArray(std::vector<float>& arr, uint Nphi, float Lb, float LX, float LY, float w_lb, float EnhanceFact);
double erfc_(double x);
float abs_(float x);
uint absInt(int n);
float factor_(unsigned int N);
float eps_scr_FormFact(float q, float ratioE, uint _NLL, float w_lb);//analitic function with FormFactor

void calc_eps_s_full_FormFact(std::vector<float>& arr, uint Nphi, float Lb, float LX, float LY, float ECl_ECr_ratio, uint _NLL, float w_lb);
//////////////////////////////////////////////////////////
//////F U N C T I O N S   D E C L A R A T I O N S /////////////////////////////
//////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
	// Here a part to send some data to python

	const uint NLL = 2;
	const uint NS = 11;
	const uint NU = 1;

	uint NE = NU * NS;
	float SZmax = 1-NE*0.5f;
	float SZmin = -0.5f*NE;

	
	//ULLONG_MAX

	uint QX_max;
	uint QY_max = 0;

	uint qqy;
	const float n2Dfrom = 1.0e11f;
	const float n2Dto = 1.9e11f;
	const float n2Dstep = 1e11f;
	//
	//Vectors of VECTORS
	std::vector<std::vector<std::complex<float>>> dataSSS;
	std::vector<std::vector<int>> indicesSSS;  //row number in each column
	std::vector<std::vector<int>> indptrSSS;
	std::vector<int> shapeSSS(2);

	uint Q_cutoff1 = 5;//cutoff of any Coulomb integral. (inclusive border)
	
	
	//-----------------------

	std::string mystringtopy = "";
	ofstream fout("file_names.txt"); //file for Boris
	fout.precision(3);
	fout.flags(ios::scientific);
	uint ijj = 0;
	int startTime, endTime, totalTime;

	startTime = (int)time(NULL);

	uint tmpQY;
	tmpQY = QY_max;

	for (float n2D = n2Dfrom; n2D <= n2Dto + 1e9; n2D += n2Dstep)
	{
		
		float LXLYratio = 1.0f;// the ratio between LY and LX
		float Lb = sqrt(float(0.5*NU / PI / n2D));
		float S = NU * NS / n2D;
		float LX = float(Lb*LXLYratio*sqrt(2 * PI*NS));
		float LY = S / LX;
		float BB = float(hbar*cc / ee / Lb / Lb);

		float ECr = float(hbar*ee*BB / cc / mm / meV_erg);
		float Ez = float(gfact*muBohr*BB / meV_erg); //Here Ez- is the full splitting between spin-states
		float ECl = float(ee*ee / eps0 / Lb / meV_erg);
		float Rydberg = float(0.5f*mm*pow(ee, 4) / hbar / hbar / eps0 / eps0 / meV_erg);
		//-----
		float width_wf = 5e-7f; //5nm - for n=1e11 - width of electron's wave function
		//width_wf = (4.98f - 0.85f*n2D / 1e11f + 0.07f*n2D*n2D / 1e22f)*1e-7f;
		float CoulombEnhanceFactor = 1.f;

		float ww_lb = width_wf / Lb;

		complex<float> zero_cmplx(0.f, 0.f);//zero complex
									
											//-----------------------------START CALCULATION--------------
		
		std::vector<float> eps_s;
		calc_eps_s_full_FormFact(eps_s, NS, Lb, LX, LY, float(ECl / ECr), NLL, ww_lb);  // here is the static dielectric function (responsible for screening). It is calculated just once and forever

		std::vector<float> FormFact_aaarr;
		calc_FormFactArray(FormFact_aaarr, NS, Lb, LX, LY, ww_lb, CoulombEnhanceFactor);//CoulombEnhanceFactor - is a compensation of Coulomb scale

		//-----  COUT All the physical parameters
		std::cout << "The list of the physical parameters:" << endl << endl;
		std::cout << "NS=" << NS << '\t' << "Nu=" << NU << '\t' << "N_LL=" << NLL << '\t' << endl << endl;
		std::cout << "n2D=" << n2D << " cm^{-2}" << '\t' << "B=" << BB * 0.0001 << "Tesla" << '\t' << "Lb=" << Lb * 1e8f << "Angstr" << endl << endl;
		std::cout << "Lx=" << LX << " cm" << '\t' << "Ly=" << LY << " cm" << '\t' << "S=" << S << "cm^2" << endl << endl;
		std::cout << "ECr=" << ECr << " meV" << '\t' << "Ez=" << Ez << " meV" << '\t' << "ECoulomb=" << ECl << " meV" << '\t' << "Rydberg=" << Rydberg << " meV" << endl << endl;



		for (float SZ = SZmax; SZ >= SZmin; SZ--)
		{
			
			uint BASIS_SIZE_Sz0;
			BASIS_SIZE_Sz0 = Basis_size_any_NU_Sz(NE, NS, NLL, SZ);

			std::vector<ullong> bas_Sz0;
			Construct_Fermi_Basis_L(bas_Sz0, NE, NS, NLL, SZ);

				

			///---------------------
			dataSSS.clear();
			dataSSS.resize(NS);
			indicesSSS.clear();
			indicesSSS.resize(NS);
			indptrSSS.clear();
			indptrSSS.resize(NS);

			/*if (SZ == -NE*0.5f)
				QY_max = 0;
			else
				QY_max = tmpQY;
			*/	QY_max = tmpQY;
			//------CHECK
			



			for (qqy = 0; qqy <= QY_max; qqy++)
			{
				QX_max = qqy + 1;
				
				shapeSSS[0] = BASIS_SIZE_Sz0;
				shapeSSS[1] = BASIS_SIZE_Sz0;

				//---------------
				std::vector<uint> blockBasisSz0(BASIS_SIZE_Sz0, 0u);

				//BE AWARE
				uint Ksi_max = (BASIS_SIZE_Sz0 > NS) ? (2 * BASIS_SIZE_Sz0) / NS : 1;//BE AWARE

				std::vector<uint> BlSizesSz0(Ksi_max, 0u);
				std::vector<int> block_signs(BASIS_SIZE_Sz0);//signs of the translated vectors



				uint Nksi0 = Filter_Basis_qy(bas_Sz0, blockBasisSz0, block_signs, BlSizesSz0, NE, NU, NS, NLL, SZ, qqy);
				//cout << "Number of classes for qy= " << qqy << " is " << Nksi0 << '\n';
				uint tmpindex;

				std::vector<uint> phi0_arr(Nksi0);//the array of indices for phi0(ksi) vectors, which are starting vectors for each class
				uint cnt = 0;

				for (uint ksi = 0; ksi < Nksi0; ksi++)
				{
					phi0_arr[ksi] = cnt;
					//CHECK momentum!
					cnt += BlSizesSz0[ksi];
				}//filling of this phi0_arr
				 //----------------------------------------------------------------

				 //---------------Matrix KSI-ETA-------------------------------------------------
				 //------------------------------------------------------------------------------

				uint TotalNonZero_elements;

				std::vector<std::complex<float>> tmpColMatr(BASIS_SIZE_Sz0); //This 1-D array is for temporary aggregation of Matrix elements in a Column
				
				int ssign = 0; // sign of the action of 4 adjacent Fermi-operators  c+c+cc	

				ullong tmpvec = 0;
				uint ind_row = 0;
				complex<float> tmp_cmplx;
				TotalNonZero_elements = 0;

				std::vector<uint> DIMqx(QX_max);//dimensions of all qx-blocks

				std::vector<uint> startIndexQxBlock(QX_max); //the start indices of qx-blocks

				std::vector<uint> nonzeroElemQx(QX_max);//the number of non-zero elements in qx-block

				uint DIM_arr = 0; //the total dimension of the full 1-D array of H_eta_ksi  
				std::vector<uint> TotalNonZeroQx(QX_max);

				std::vector<uint> ksi_counter(QX_max);
				uint eta_counter = 0;

				for (uint qx = 0; qx < QX_max; qx++) //here we just calculate the dimensions of the qx-blocks for matrices H_eta_ksi(qx)
				{
					TotalNonZeroQx[qx] = 0;//set it zero-valued array
					ksi_counter[qx] = 0;
					DIMqx[qx] = 0;//start from zero
					startIndexQxBlock[qx] = DIM_arr;//the start indices of successive qx-blocks
					for (uint i = 0; i < Nksi0; i++)
					{
						if (((qx * BlSizesSz0[i]) % NS) == 0)
						{
							DIMqx[qx]++;//this variable stores the dimension of the ksi-eta matrix for definite qx-value 
							DIM_arr++;
						}
					}
				}//-----------------------------------------------------

				 //This joint 1-D array is for temporary aggregation of Matrix elements numerated by "eta" and for all qx-projections
				std::vector<std::complex<float>> tmpEtaElements(DIM_arr);//full length joined for all qx-blocks
																		 //CHECK that the array is zero-valued

				uint pos;  //the position in the Block-qxqy basis
				uint phi0_ksi;//index of the phi0 vector for each class "ksi"
				uint Bsize_ksi;
				uint Bsize_eta;
				complex<float> tSum;

				for (uint qx = 0; qx < QX_max; qx++)
				{
					indptrSSS[qx].push_back(0);   /////SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
				}


				uint ksi_ind_original;
				uint eta_ind_original;
				
				uint eta_indexBlock;
				uint ksi_indexBlock;
				


				for (uint ksi = 0; ksi < Nksi0; ksi++)//ksi_i - index of class
				{
					//first set diagonal element:
					std::cout << "Sz=" << SZ << "qy=" << qqy << "  => ksi = " << ksi << '\n';
					phi0_ksi = phi0_arr[ksi];
					ksi_ind_original = blockBasisSz0[phi0_ksi];

					tmpColMatr[ksi_ind_original] = CalcSingleParticle_En(bas_Sz0[ksi_ind_original], NLL, NS, ECr, Ez); //col=row for diagonal element



					for (uint i1 = 0; i1 < NS; i1++)
						for (uint i2 = 0; i2 < NS; i2++)
							for (uint i3 = 0; i3 < NS; i3++)
								for (uint i4 = 0; i4 < NS; i4++)
									for (uint n1 = 0; n1 < NLL; n1++)
										for (uint n2 = 0; n2 < NLL; n2++)
											for (uint n3 = 0; n3 < NLL; n3++)
												for (uint n4 = 0; n4 < NLL; n4++)
													for (uint s14 = 0; s14 < 2; s14++)
														for (uint s23 = 0; s23 < 2; s23++)
														{
															assert(bas_Sz0.size() > ksi_ind_original);

															if ((tmpvec = TwoParticleMixingOperator(bas_Sz0[ksi_ind_original], ssign, i1, i2, i3, i4, n1, n2, n3, n4, s14, s23, NS)) > 0) //it means that action of the operator on the vector is non-zero
															{
																ind_row = index_vector(tmpvec, bas_Sz0, BASIS_SIZE_Sz0);

																tmp_cmplx = float(ssign)*MatrixElement((int)i1, (int)i2, (int)i3, (int)i4, n1, n2, n3, n4, NS, Lb, LX, LY, ECl, eps_s, FormFact_aaarr, Q_cutoff1);

																assert(tmpColMatr.size() > ind_row);
																tmpColMatr[ind_row] += tmp_cmplx;
															}
														}

					//Here we have the column "ksi" of the matrix ready for calculating the matrix elements H_ksi_eta(qx) for the set of momenta "qx" in array qxarr

					for (uint row = 0; row < BASIS_SIZE_Sz0; row++)
						round_cmplx(tmpColMatr[row]);

					//////////////////////////////////////////////////
					//Calculating and recording the matrix elements H_eta_ksi to multiple qx-files

					//------------------------
					Bsize_ksi = BlSizesSz0[ksi]; //Block_size of "ksi"
												 //for (uint qx = 0; qx < QX_max; qx++)
												 //{
					for (uint qx = 0; qx < QX_max; qx++)
					{
						if (((qx*Bsize_ksi) % NS) == 0) //the condition of this ksi-block to participate
						{//OK here the qx-test is passed
							int numberss = 0;///SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSs

							eta_counter = 0;

							for (uint ieta = 0; ieta < Nksi0; ieta++)//filling the matrix elements H_ksi_eta(qx) for the set of momenta "qx" in array qxarr
							{
								Bsize_eta = BlSizesSz0[ieta]; //Block_size of "ieta"-block

								if (((qx*Bsize_eta) % NS) == 0)//it means that this "qx" is compatible with both "ksi" and "eta"-block indices
								{//so now we calculate H_eta_ksi

									tSum = zero_cmplx;

									for (uint k = 0; k < Bsize_eta; k++)
									{
										eta_indexBlock = phi0_arr[ieta] + k;   //this is the successive index of our element "phi_eta(j)" in the array "blockBasisSz0"
										eta_ind_original = blockBasisSz0[eta_indexBlock];//here we get an index of the vector in the original array "Bas_Sz0" 

										tmp_cmplx = tmpColMatr[eta_ind_original];
										
										ssign = block_signs[eta_indexBlock];

										tSum = tSum + tmp_cmplx * polar<float>(1.0f, float(-2.0f*PI*qx*k / NS))*float(ssign);
									}

									tmp_cmplx = sqrt(1.0f* Bsize_ksi / Bsize_eta)*tSum;//here we get the true matrix element H_eta_ksi and are ready to save it to the output matrix-file
									round_cmplx(tmp_cmplx);

									if (tmp_cmplx != zero_cmplx)
									{
										TotalNonZero_elements++;
										TotalNonZeroQx[qx]++;

										numberss++;
										dataSSS[qx].push_back(tmp_cmplx);
										indicesSSS[qx].push_back(eta_counter);/////////////////////////////////////////////
																			  //WRITE ALSO TO THE MATRIX FILE							
									}
									eta_counter++;//increment the actual eta counter
								}
							}//for "eta" - end
							indptrSSS[qx].push_back(indptrSSS[qx].back() + numberss);//SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS						
						}//if  qx-check
						 //otherwise this qx-value is not compatible with the current "ksi"
						 // 
						ksi_counter[qx]++;
						//cout << "qy,qx= (" << qqy << "," << qx << ")	eta=" << eta_counter << "ksi= " << ksi_counter[qx] << '\n';
					}//end of qx-cycle



					 //ERASING the array "tmpColMatr"
					for (uint row = 0; row < BASIS_SIZE_Sz0; row++)
					{
						tmpColMatr[row] = 0;//erase all matrix elements in the array for the next cycle iteration
					}
				}

				// END OF Sz=0 MATRIX  ksi and eta MATRIX
				for (uint qx = 0; qx < QX_max; qx++)
				{
					//std::cout << "Dim_ksi = " << ksi_counter[qx] << " , DIM_eta=" << eta_counter << '\n';
					//std::cout << "Total Non-zero= " << TotalNonZeroQx[qx] << '\n';

					ijj++;//file number

					std::string copyofstring;
					if (SZ == 0) copyofstring = "00";
					else
						copyofstring = std::to_string(int(round(SZ)));

					std::string iternumb;
					if (ijj < 10)
						iternumb = "00" + std::to_string(ijj);
					else
						if ((ijj >= 10) && (ijj <100))
							iternumb = "0" + std::to_string(ijj);
						else
							iternumb = std::to_string(ijj);



					mystringtopy.clear();
					mystringtopy = iternumb + "_qx_" + std::to_string(qx) + "_qy_" + std::to_string(qqy) + "_s_" + copyofstring + "_NS_" + std::to_string(NS) + ".npz";
					fout << mystringtopy + "   " + std::to_string(n2D) << endl;
					shapeSSS[0] = DIMqx[qx];////SSSSS  
					shapeSSS[1] = DIMqx[qx];////SSSSS 
					cnpy::npz_save(mystringtopy, "shape", &shapeSSS[0], { 2 }, "w"); //"a" appends to the file we created above
					cnpy::npz_save(mystringtopy, "data", &dataSSS[qx][0], { TotalNonZeroQx[qx] }, "a"); //"a" appends to the file we created above
					cnpy::npz_save(mystringtopy, "indices", &indicesSSS[qx][0], { TotalNonZeroQx[qx] }, "a"); //"a" appends to the file we created above
					cnpy::npz_save(mystringtopy, "indptr", &indptrSSS[qx][0], { indptrSSS[qx].size() }, "a"); //"a" appends to the file we created above

					dataSSS[qx].clear();
					indicesSSS[qx].clear();
					indptrSSS[qx].clear();


				}
				
				//cout << "everything is clear" << endl;
				std::cout << "FINISHED  (qy)= " << qqy << endl;


			}//end of qqy-cycle
						
		}//end of Sz-cycle

		cout << "n2D= " << n2D << "  finished. Ecr-Ezm = " << ECr - Ez << "meV\n";
	}
	//-----------
	endTime = (int)time(NULL);
	totalTime = endTime - startTime;

	std::cout << "Total time = " << totalTime << "seconds" << '\n';
	
	fout.close();
	//M_file.close();
	int ssign = 0;
	std::cin >> ssign;

	return 0;
	//-----------------

}

/////////////////////////////////////////////
/////////////////////////////////////////////
/////////////////////////////////////////////
////F U N C T I O N S////////////////////////
/////////////////////////////////////////////
////F U N C T I O N S////////////////////////
/////////////////////////////////////////////
////F U N C T I O N S////////////////////////




void Construct_Fermi_Basis_L(std::vector<ullong>& s_Sz, uint Ne, uint NS, uint NLL, float Sz)
{
	/// the constructor of the basis states of electrons
	
	uint Sz_basis_DIM = Basis_size_any_NU_Sz(Ne, NS, NLL, Sz);
	//ullong* s_Sz = new ullong[Sz_basis_DIM];
	s_Sz.resize(Sz_basis_DIM);//here must be +++++++++++++++++1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
	s_Sz[0] = 0;
	ullong tmp = 0x1;
	//
	uint counter_Sz = 0;//set counter to 0. This is to become a dimension of Basis-array
						/// ---
	for (uint i = 0; i<Ne - 1; i++)
	{
		tmp = tmp << 1;
		tmp++;
	}

	ullong vect = tmp;  //stores a previous vector for calculation of a new one
						//the first basis state
	if (Sz_full_represent(tmp, NS, NLL) == Sz)  //check spin and choose state
	{
		s_Sz[counter_Sz] = tmp;				//check spin and choose state
		counter_Sz++;
	}

	ullong bbit;

	while (counter_Sz < Sz_basis_DIM)
	{
		bbit = ullong(0x1);
		int i = 0;
		while ((tmp = (vect&bbit)) != 0)
		{
			bbit = bbit << 1;
			i++;
		}

		if (i >= 1)
		{
			vect = vect + (ullong(0x1) << (i - 1));
		}
		else// search for the next "0"-symbol
		{
			int k1 = 0;
			while ((tmp = (vect&bbit)) == 0)
			{
				bbit = bbit << 1;
				k1++;
			}
			int k2 = k1;
			while ((tmp = (vect&bbit)) != 0)
			{
				bbit = bbit << 1;
				k2++;
			}
			vect = vect + (ullong(0x1) << k1) + (ullong(0x1) << (k2 - k1 - 1)) - 1;
		}

		if (Sz_full_represent(vect, NS, NLL) == Sz)  //check spin and choose state
		{
			assert(s_Sz.size() > counter_Sz);
			
			s_Sz[counter_Sz] = vect;				//check spin and choose state
			counter_Sz++;
		}
	}
	//-----------		
}
//-----------												


//-----------

uint Combinat_N_m(uint N, uint m)
{
	if (m > N)
		return 0;
	if ((m == N) || (m==0))
		return 1;
	//else
	float tmp = 1.0f;
	for (int j = 0; j < m; j++)
	{
		tmp = tmp*(N - j) / (m - j);	
	}
	return (uint)round(tmp);
}



//-----------
uint Basis_size_any_NU_Sz(uint _Ne, uint _NS, uint _NLL, float Sz) 
{
	uint n1 = uint(_Ne*0.5f + Sz);
	uint n2 = uint(_Ne*0.5f - Sz);
	

	return Combinat_N_m(_NLL*_NS, n1)*Combinat_N_m(_NLL*_NS, n2);
}

//-----------

//OLD BASIS DIM   FORMULA

float CalcSingleParticle_En(ullong vector, uint _NLL, uint _NS, float Ecr, float Ezm)
{
	ullong tmp = 0x1;
	float res = 0;

	for (uint n = 0; n<_NLL; n++)
	{
		for (uint i = 0; i<_NS; i++) //Spin-down counter is the first for positive g-factor
		{
			if (vector&tmp)
				res = res + (n + 0.5f)*Ecr - Ezm * 0.5f;
			tmp = tmp << 1;
		}
		for (uint i = 0; i<_NS; i++) // Spin-up counter is here
		{
			if (vector&tmp)
				res = res + (n + 0.5f)*Ecr + Ezm * 0.5f;
			tmp = tmp << 1;
		}
	}
	return res;
	//-------
}

//-----------old    U INT------------------------

///////////////////////////
void display_vect(ullong vv, uint NS, uint NLL)
{
	ullong bbit = 0x1;
	ullong tmp;
	uint i = 0; //just counter

	for (uint n = 0; n < NLL; n++)
	{
		std::cout << "/";
		for (uint Sz = 0; Sz <= 1; Sz++)
		{
			std::cout << "-";
			for (uint p = 0; p < NS; p++)
			{
				tmp = (vv&bbit) >> i;
				i++;
				bbit = bbit << 1;
				std::cout << tmp;
			}
			std::cout << "-";
		}
		std::cout << "/";
	}

	std::cout << '\n';
}
///////////////////////////
///////////////////////////
ullong Birth(ullong vect, uint pos)
{

	if ((vect&(ullong(0x1) << pos)) != 0)
		return 0;  //here Birth is failed.  return 0
	else
	{
		return vect + (ullong(0x1) << pos);
	}
}
///////////////////////////
///////////////////////////
ullong Death(ullong vect, uint pos)
{

	if ((vect&(ullong(0x1) << pos)) != 0)
		return vect - (ullong(0x1) << pos);
	else
		return 0;
}
///////////////////////////

int SignBirth(ullong vector, uint pos)
{
	ullong bbit = 0x1;
	uint sum = 0;
	for (uint i = 0; i<pos; i++)
	{
		if (vector & bbit)
			sum++;
		bbit = bbit << 1;
	}
	if (sum & 0x1) //i.e. sum is odd
		return -1;
	else
		return 1;
}
/////////
//////////////////////
int SignDeath(ullong vector, uint pos)
{
	ullong bbit = 0x1;
	uint sum = 0;
	for (uint i = 0; i<pos; i++)
	{
		if (vector & bbit)
			sum++;
		bbit = bbit << 1;
	}
	if (sum & 0x1) //i.e. sum is odd
		return -1;
	else
		return 1;
}
//////////////////////
ullong TwoParticleMixingOperator(ullong vector, int &ssign, uint i1, uint i2, uint i3, uint i4, uint n1, uint n2, uint n3, uint n4, uint spin14, uint spin23, uint Nphi)
{
	ullong result = vector;
	uint pos1 = Nphi * (2 * n1 + spin14) + i1;
	uint pos2 = Nphi * (2 * n2 + spin23) + i2;
	uint pos3 = Nphi * (2 * n3 + spin23) + i3;
	uint pos4 = Nphi * (2 * n4 + spin14) + i4;
	int sign = 1;  // Sign of the Fermi-operators

	if ((result = Death(result, pos4)) != 0)
	{
		sign = sign * SignDeath(result, pos4);
		//display_vect(result,NS,NLL);//---
		if ((result = Death(result, pos3)) != 0)
		{
			sign = sign * SignDeath(result, pos3);
			//display_vect(result,NS,NLL);//--
			if ((result = Birth(result, pos2)) != 0)
			{
				//display_vect(result,NS,NLL);//-
				sign = sign * SignDeath(result, pos2);

				if ((result = Birth(result, pos1)) != 0)	 //here is the ultimate vector
				{
					sign = sign * SignDeath(result, pos1);		// here is the ultimate "sign" which we use in the end of calculation V(...)
																//display_vect(result,NS,NLL); 
																//--OK, we have a non-zero answer
					ssign = sign;
					return result;
					//-
				}
			}
		}
	}
	ssign = 0;//for zero-answer
	return result; //here is the ultimate "zero-answer
}


//////////////////////
uint DeltaSymbol(int pos1, int pos2, uint Nphi)
{
	int p1 = (pos1 + Nphi) % Nphi;
	int p2 = (pos2 + Nphi) % Nphi;

	if (p1 != p2)
		return 0;
	return 1;
}


//--------
float sq(float x)
{
	return x * x;
}
//-----------
int sign_u(uint k, uint m)
{
	if (k<m) return -1;
	if (k>m) return 1;
	return 0;//otherwise
}
//-----------
uint factorial(uint N)
{
	if (N <= 1)
		return 1;
	else return N * factorial(N - 1);
}
//---------


//---------
float FFnmReal(float x, uint n, uint m)
{
	int n1 = (int)n;
	int m1 = (int)m;

	if (n1<m1)
	{
		int tmp = m1;
		m1 = n1;
		n1 = tmp;
	}          //So  here  we set n1>=m1

	float Laguerr;
	int n_m = n1 - m1;


	switch (m1)
	{
	case 0:	Laguerr = 1; break;
	case 1: Laguerr = 1 + n_m - x; break;
	case 2: Laguerr = (n_m + 2)*(n_m + 1)*0.5f - (n_m + 2)*x + 0.5f*x*x; break;
	default: Laguerr = 0; break; //We avoid calculating at LL>=3 occupied
	}

	return sqrt(1.0f*factorial(uint(m1)) / factorial(uint(n1)))*exp(-x * 0.5f)*Laguerr;
}
//----------------
//This function applies for NON-ZERO argument ONLY
uint index_vector(ullong vect, const std::vector<ullong>& arr_vect, uint NMAX)
{
	uint left = 0;
	uint right = NMAX - 1;
	uint ind = (left + right) / 2;

	ullong tmp;
	if (vect < 0xffffffffffffffff)
	{
		do
		{
			tmp = arr_vect[ind];

			if (tmp < vect)
			{
				left = ind;
				ind = (left + right + 1) / 2;  // +1 - to avoid edge-effect
			}
			if (tmp > vect)
			{
				right = ind;
				ind = (left + right) / 2;
			}

		} while (tmp != vect);
		//MISTAKE!!!
		return ind;
	}
	else return 0;
}
////-----------------------
void round_cmplx(complex<float> &z)
{
	complex <float> tmp(z);
	if (abs(z.real())<float_err)
	{
		tmp.real(0);
	}
	if (abs(z.imag())<float_err)
	{
		tmp.imag(0);
	}
	z = tmp;
	return;
}
//---////////////////


complex<float> MatrixElement(int i1, int i2, int i3, int i4,
	uint n1, uint n2, uint n3, uint n4,
	uint Nphi, float Lb, float LX, float LY, float ECl, const std::vector<float>& eps_s, const std::vector<float>& FormFactArray, uint Q_cutoff)
{
	//// Now calculate the Coulomb factor V(i1,i2,i3,i4), where an integration over qx,qy is performed
	complex<float> SumV;
	float q2, qlx, qly;
	//float qx_factor=(float)(2*PI*Lb/LX*(i3-i1));
	float tt = 0;

	complex<float> qxy_n1n4(0.f, 0.f);
	complex<float> qxy_n2n3(0.f, 0.f);



	complex<float> V_MatrElem(0.f, 0.f); //Sum-variable
	int Q_cut = (int)Q_cutoff;
	int qLim = (int)Nphi;
	uint indexInArray;

	for (int qx = -qLim + 1; qx<qLim; qx++)
	{
		SumV = complex<float>(0.f, 0.f);

		for (int qy = -qLim + 1; qy<qLim; qy++)
		{
			if ((DeltaSymbol(i1, i4 + qy, Nphi)*DeltaSymbol(i2, i3 - qy, Nphi)) && ((qx != 0) || (qy != 0)))  //First check Cronecker symbols and that q!=0
			{
				if ((absInt(qx) <= Q_cut) && (absInt(qy) <= Q_cut))
				{
					qlx = float(2 * PI*Lb)*qx / LX;
					qly = float(2 * PI*Lb)*qy / LY;
					q2 = sq(qlx) + sq(qly);

					//-----------------
					//-----------polynom of qx and qy
					if (n1 != n4)  //otherwise it equals 1
					{
						qxy_n1n4 = complex<float>(float(sign_u(n1, n4)*qly / sqrt2), float(qlx / sqrt2));
						qxy_n1n4 = pow_cmplxint(qxy_n1n4, abs((int)n1 - (int)n4));
					}
					else
						qxy_n1n4 = complex<float>(1.0f, 0.f);

					if (n2 != n3)  //otherwise it equals 1
					{
						qxy_n2n3 = complex<float>(float(-sign_u(n2, n3)*qly / sqrt2), float(-qlx / sqrt2));
						qxy_n2n3 = pow_cmplxint(qxy_n2n3, abs((int)n2 - (int)n3));
					}
					else
						qxy_n2n3 = complex<float>(1.0f, 0.f);

					indexInArray = Nphi * absInt(qy) + absInt(qx);
					SumV = SumV + FormFactArray[indexInArray] / sqrt(q2) / eps_s[indexInArray] * FFnmReal(0.5f*q2, n1, n4)*FFnmReal(0.5f*q2, n2, n3)*qxy_n1n4*qxy_n2n3;
				}
			}
		}
		V_MatrElem = V_MatrElem + SumV * polar<float>(1.0f, float(2.0f*PI*qx*(i3 - i1) / Nphi));
	}
	return 0.5f / Nphi * ECl*V_MatrElem;  //where EC-interparticle Coulomb energy;
}
//--------------




//----------------------------
//----------------------------
float Sz_full_represent(ullong vector, uint _Nphi, uint _NLL) // it means that both spin states of each LL are considered
																  //APPLICABLE ONLY FOR EVEN NUMBER OF ELECTRONS!!!
{
	int S_up = 0, S_down = 0;
	ullong tmp = 0x1;

	for (uint n = 0; n<_NLL; n++)
	{
		for (uint i = 0; i<_Nphi; i++) //Spin-down counter is the first for positive g-factor
		{
			if (vector&tmp)
				S_down++;
			tmp = tmp << 1;
		}
		for (uint i = 0; i<_Nphi; i++) // Spin-up counter is here
		{
			if (vector&tmp)
				S_up++;
			tmp = tmp << 1;
		}
	}

	return 0.5f*(S_up - S_down);
}
//--------------------------------

complex<float> pow_cmplxint(complex<float> z, uint t)
{
	complex<float> unity(1.0f, 0.0f);
	if (t == 0)
		return unity;
	if (t == 1)
		return z;
	else
	{
		return polar<float>(pow(abs<float>(z), (int)t), t*arg<float>(z));
	}
}

//---------------------------------------------------------------------

//-----------------
complex<float> conjugate(complex<float> z)
{
	return complex<float>(z.real(), -z.imag());
}
//-----------------
uint absInt(int n)
{
	if (n < 0) return uint(-n);
	else
		return uint(n);
}
//-----------------


ullong TranslVectorSigned(ullong vector, int &ssign, uint NS, uint NLL)//here we perform a collinear translation of each spin sub-level by 1 node
{
	const ullong ulx1 = 0x1;
	ullong result = 0;
	ullong tmp;
	uint pos0 = 0;
//	uint posMax;
	//ullong MASK;
	ullong bbit = ulx1;
	ullong blockNS = ulx1;



	for (int j = 0; j < int(NS) - 1; j++)
		blockNS = (blockNS << 1) + 1;

	int sign = 1;//original sign
	ullong tmpvect = vector;
	uint newpos;
	
	for (uint i = 0; i<NLL; i++)
		for (uint Sz = 0; Sz <= 1; Sz++)
		{
			pos0 = (2 * i + Sz)*NS;
			bbit = ulx1 << pos0;
			for (uint j = 0; j < NS; j++)
			{
				if (vector&bbit)
				{
					tmpvect = Death(tmpvect, pos0 + j);
					sign = sign * SignDeath(tmpvect, pos0 + j);
					
				}
				bbit = bbit << 1;//shift
			}

		}
	
	for (int i = NLL - 1; i >= 0; i--)
		for (int Sz = 1; Sz >= 0; Sz--)
		{
			//Create with shift
			pos0 = (2 * uint(i) + uint(Sz))*NS;
			bbit = ulx1 << (int(pos0 + NS) - 1);       //the highest node on the last spin-sublevel
			for (int j = int(NS) - 1; j >= 0; j--)
			{
				if (vector&bbit)
				{
					newpos = pos0 + (uint(j) + 1) % NS;
					tmpvect = Birth(tmpvect, newpos);
					sign = sign * SignBirth(tmpvect, newpos);					
				}
				bbit = bbit >> 1;//shift back
			}
		}

	ssign = sign;
	return tmpvect;
}
//------------------------------

//------------------------------
uint Calcul_qy_for_vector(ullong vector, uint Ne, uint NS, uint NLL)
{
	uint sum = 0; // for the overall answer
	uint pos;
	ullong tmp;
	uint count;

	for (uint i = 0; i < NS; i++)
	{
		count = 0;
		for (uint n = 0; n<NLL; n++)
			for (uint sz = 0; sz <= 1; sz++)
			{
				pos = (2 * n + sz)*NS + i;  //position of the electron in the vector
				tmp = ullong(0x1) << pos;
				if (tmp & vector)
					count++;
			}     //end of counter of the electrons in the i-th node
		sum += count * i;
	}
	sum = sum + (Ne*(NS + 1)) / 2;
	return sum % NS; //mod NS. So the answer will always be in the range 0..(NS-1)

}

//------------------------------
//uint Filter_Basis_qy(const std::vector<uint>& in_basSz, std::vector<uint>& out_blockBasis, std::vector<int>& out_signs, std::vector<uint>& out_sizes,
//	uint Ne, uint NU, uint NS, uint NLL, float Sz, uint qy); //create sub_basis of states with the fixed qy-projection
//------------------------------
uint Filter_Basis_qy(const std::vector<ullong>& in_basSz, std::vector<uint>& out_blockBasis, std::vector<int>& out_signs, std::vector<uint>& out_sizes, uint Ne, uint NU, uint NS, uint NLL, float Sz, uint qy) //create sub_basis of states with the fixed Sz projection
{//returns number of CLASSES, satisfying the input conditions

 /// the constructor of the basis states of electrons with filters Sz=Sz,  qy=qy
	uint Sz_basis_DIM;
	Sz_basis_DIM = Basis_size_any_NU_Sz(Ne, NS, NLL, Sz);
	
	std::vector<ullong> tmp_in(Sz_basis_DIM);
	
	if ((qy < 0) || (qy >= NS))
		return 0; //ERROR



	for (int j = 0; j < Sz_basis_DIM; j++)
	{
		tmp_in[j] = in_basSz[j];
	}//temporary implementation of the input basis
	 //////////////////////
	uint ksi = 0;
	uint counter = 0;
	ullong phiStart;   //proizvodyaschaya funktziya
	ullong ttt;
	uint ind_tmp; //temporary variable for index
	int signTrans = 1;
	int sign = 1;

	for (uint i = 0; i < Sz_basis_DIM; i++)
	{
		if ((tmp_in[i]) && (Calcul_qy_for_vector(tmp_in[i], Ne, NS, NLL) == qy)) //vector is non-zero and its momentum equals qy
		{
			signTrans = 1;
			phiStart = tmp_in[i];
			out_blockBasis[counter] = i;
			out_signs[counter] = signTrans;//just for start vector
			out_sizes[ksi] = out_sizes[ksi] + 1;
			counter++;
			//cyclic Translations:
			ttt = phiStart;

			//cout << "ksi=" << ksi << endl;//--------------------------
			//cout << "\t PhiStart= ";//-----------------------------
			//display_vect(phiStart, NS, NLL);//-------------------------

			while ((ttt = TranslVectorSigned(ttt, sign, NS, NLL)) != phiStart)
			{
				//display_vect(ttt, NS, 4);
				//std::cout.write(reinterpret_cast<const char*>(&ttt), sizeof ttt);
				//std::cout << ttt <<endl;//------------------------

				ind_tmp = index_vector(ttt, in_basSz, Sz_basis_DIM);

				signTrans = signTrans * sign;
				out_blockBasis[counter] = ind_tmp;
				out_signs[counter] = signTrans;
				out_sizes[ksi] = out_sizes[ksi] + 1;
				tmp_in[ind_tmp] = 0; //make it zero to clean the input array
				counter++;
			}//so here the cyclic group of class "ksi"  is finished			
			ksi++;
		}
	}
	
	return ksi;//Everything is OK. We return number of classes
}
//-----------END of Construct_Fermi_basis_SPIN



float abs_(float x)
{
	if (x >= 0) return x;
	else return -x;
}
//////////////////////////////
float factor_(unsigned int N)
{
	if (N<2)  return 1.f;
	else return N * factor_(N - 1);
}
//---------------------------------


float FormFact(float q, float w_lb)
{
	if (q == 0)
		return 1.f;
	else
		return exp(sq(q*w_lb))*(float)erfc_(double(q*w_lb));
}

//--------------

void calc_FormFactArray(std::vector<float>& arr, uint Nphi, float Lb, float LX, float LY, float w_lb, float EnhanceFact)
{
	arr.resize(Nphi*Nphi);

	for (uint qx = 0; qx<Nphi; qx++)
		for (uint qy = 0; qy<Nphi; qy++)
		{
			arr[qy*Nphi + qx] = EnhanceFact * FormFact(float(2 * PI*Lb*sqrt(sq(qx / LX) + sq(qy / LY))), w_lb);
		}	
}

///-/-------------------------
float eps_scr_FormFact(float q, float ratioE, uint _NLL, float w_lb)
{
	float sum = 0;

	if (q != 0)
	{
		for (int n = _NLL; n <= 9; n++)
		{
			sum = sum + pow(q, (int)(2 * n - 1)) / n / factor_(n)*pow(0.5, (int)n);
		}
		return 1 + 2* ratioE*exp(-0.5*q*q)*sum*FormFact(q, w_lb);//for NU=1
	}
	else
	{
		return 1.f;
	}
}
//---------------------------------

//--------------
void calc_eps_s_full_FormFact(std::vector<float>& arr, uint Nphi, float Lb, float LX, float LY, float ECl_ECr_ratio, uint _NLL, float w_lb)
{
	arr.resize(Nphi*Nphi);

	for (uint qx = 0; qx<Nphi; qx++)
		for (uint qy = 0; qy<Nphi; qy++)
		{
			arr[qy*Nphi + qx] = eps_scr_FormFact(float(2 * PI*Lb*sqrt(sq(qx / LX) + sq(qy / LY))), ECl_ECr_ratio, _NLL, w_lb);
		}	
}
//------------------
//--------------
float* calc_FormFactArray(uint Nphi, float Lb, float LX, float LY, float w_lb)
{
	float* arr;

	arr = new float[Nphi*Nphi];

	for (uint qx = 0; qx<Nphi; qx++)
		for (uint qy = 0; qy<Nphi; qy++)
		{
			arr[qy*Nphi + qx] = FormFact(float(2 * PI*Lb*sqrt(sq(qx / LX) + sq(qy / LY))), w_lb);
		}
	return arr;
}
//////////////////////////////
//------
uint min(uint a, uint b)
{
	if (a < b)
		return a;
	else return b;
}
//-----------SINGLE-PARTICLE diagonal
uint max(uint a, uint b)
{
	if (a > b)
		return a;
	else return b;
}
///////////////////////////////////////////
double erfc_(double x)
{
	double sum = 0;
	const int N = 18;

	if (x <= 3)
	{
		for (int k = 0; k<N; k++)
		{
			sum = sum +
				pow(x, (int)(4 * k + 1)) / factor_(2 * k) / (4 * k + 1)
				- pow(x, (int)(4 * k + 3)) / factor_(2 * k + 1) / (4 * k + 3);
		}

		return 1 - sum * 1.12837917;
	}
	else  //asymptotical form
	{
		for (int k = 0; k<4; k++)
		{
			sum = sum + factor_(4 * k) / factor_(2 * k) / pow(2 * x, (int)(4 * k)) - factor_(4 * k + 2) / factor_(2 * k + 1) / pow(2 * x, 4 * k + 2);
		}
		return sum * exp(-x * x) / x / 1.77245;
	}
}