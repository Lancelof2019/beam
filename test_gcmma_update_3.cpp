#include <gcmma/GCMMASolver.h>
#include <mma/MMASolver.h>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "cc_membrane.cc"

double Squared(double x) { return x*x; }

struct Problem {
	int n, m;
	std::vector<double> x0, xmin, xmax;
	double x_density_percell;
    double *cell_arry;
	int dofs_all_num;
	double norm_value;
	double x_density_coff;
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      deallog.depth_console(0);
      Parameters::AllParameters parameters("parameters.prm");
      if (parameters.automatic_differentiation_order == 0)
        {
          std::cout << "Assembly method: Residual and linearisation are computed manually." << std::endl;

          // Allow multi-threading
          Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv,
                                                              dealii::numbers::invalid_unsigned_int);

          typedef double NumberType;
          Solid<dim,NumberType> solid_3d(parameters);
          solid_3d.run();
        }
	   
       double problem_mu=solid_3d.mu;//density of object
	   
	   int active_cells=solid_3d.dofs_per_cell;//numbers of cells
	   
	   dofs_all_num=solid_3d.tangent_matrix.m();
	   
	   Eigen::SparseMatrix<double> eigen_tangent_matrix(solid_3d.tangent_matrix.m(), solid_3d.tangent_matrix.n());
       #pragma omp parallel for
       for (unsigned int i = 0; i < solid_3d.tangent_matrix.n(); ++i) {
          for (dealii::BlockSparseMatrix<double>::const_iterator it = solid_3d.tangent_matrix.begin(i);
             it !=solid_3d.tangent_matrix.end(i); ++it) {
             eigen_tangent_matrix.insert(it->row(), it->column()) = it->value();
          }
        }
		
		
		Eigen::MatrixXd dense_matrix = eigen_tangent_matrix.toDense();
        Eigen::MatrixXd inv_dense_matrix = dense_matrix.inverse();
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(dense_matrix);
        
		if (eigensolver.info() != Eigen::Success) {
           std::cout << "Failed to compute eigenvalues." << std::endl;
        return 1;
       }
        /*Eigen::VectorXd eigenvalues = eigensolver.eigenvalues();
        Eigen::MatrixXd eigenvectors = eigensolver.eigenvectors();
		Eigen::MatrixXd B = eigenvectors;
        Eigen::MatrixXd D = Eigen::MatrixXd::Zero(B.cols(), B.cols());
		#pragma omp parallel for
        for (int i = 0; i < B.cols(); i++) {
          D(i, i) = std::sqrt(eigenvalues(i));
        }
		   Eigen::VectorXd solver_vector(solid_3d.solution_n.size());

       #pragma omp  parallel for
       for (unsigned int i = 0; i < solid_3d.solution_n.size(); ++i) {
             solver_vector(i) = solid_3d.solution_n[i];
         }
       Eigen::VectorXd solver_gradient=-1*inv_dense_matrix*D*solver_vector;
	   */
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<





	Problem()
		: n(active_cells)//number of cells or grids
		, m(2)//destination file with constraint :1)desination 2)v<=volumn_value
		, x_density_coff(int density_coff)//got density coefficient from outside:1*density_coff
		, x0({4, 3, 2})
		, xmin(n, 0.0)
		, xmax(n, 5.0)
	{ 	  
	
	   x_density_percell=x_density_coff*problem_mu;
	   cell_array=new double[n];
	   #pragma omp parallel for
	   for(int i=0;i<n;i++){
		   cell_array[i]=x_density_percell;
	   }
	
	}

	/*void Obj(const double *x, double *f0x, double *fx) {
		f0x[0] = 0;
		for (int i = 0; i < n; ++i) {
			f0x[0] += x[i]*x[i];
		}
		fx[0] = Squared(x[0] - 5) + Squared(x[1] - 2) + Squared(x[2] - 1) - 9;
		fx[1] = Squared(x[0] - 3) + Squared(x[1] - 4) + Squared(x[2] - 3) - 9;
	}
	*/
	void Obj(const double *x, double *f0x, double *fx) {
		//J(x)=norm1(x);x is slover from Ax=b;
		
		
		
		for(int i=0;i<dofs_all_num;i++){
			
			//f0x[i]=dense_matrix[i];
			fx[0]+=abs(f0x[i]);
		}
				
				
				
		for(int j=0;j<;j++){
			if(fx[1]<=0){
			fx[1]+=cell_array[i]-0.5*solid_3d.mu;
			}
		}
				
	   }


	void ObjSens(const double *x, double *f0x, double *fx, double *df0dx, double *dfdx) {
		
		Obj(x, f0x, fx);
	  
	    Eigen::VectorXd df0dx_vector(dofs_all_num);
	    for(int i=0;i<dofs_all_num;i++){
			
			df0dx[i]=f0x[i]/fx[0];
			df0dx_vector[i]=df0dx[i];
		}//求出目标函数导数df0dx
			
		
		//>>>伴随法
		Eigen::VectorXd eigenvalues = eigensolver.eigenvalues();
        Eigen::MatrixXd eigenvectors = eigensolver.eigenvectors();
		Eigen::MatrixXd B = eigenvectors;
        Eigen::MatrixXd D = Eigen::MatrixXd::Zero(B.cols(), B.cols());
		#pragma omp parallel for
        for (int i = 0; i < B.cols(); i++) {
          D(i, i) = std::sqrt(eigenvalues(i));
        }
		
	   Eigen::VectorXd solver_vector(solid_3d.solution_n.size());

       #pragma omp  parallel for
       for (unsigned int i = 0; i < solid_3d.solution_n.size(); ++i) {
             solver_vector(i) = solid_3d.solution_n[i];
         }
       Eigen::VectorXd solver_gradient=-df0dx_vector*inv_dense_matrix*D*solver_vector;
	   
	   for(int i=0;i<dofs_all_num;i++){
		   
		   dfdx[i]=solver_gradient[i];
		   
	   }//求出梯度
		
	}
	
		/*void ObjSens(const double *x, double *f0x, double *fx, double *df0dx, double *dfdx) {
		Obj(x, f0x, fx);
		for (int i = 0; i < n; ++i) {
			df0dx[i] = 2*x[i];
		}
		int k = 0;
		dfdx[k++] = 2 * (x[0] - 5); dfdx[k++] = 2 * (x[0] - 3);
		dfdx[k++] = 2 * (x[1] - 2); dfdx[k++] = 2 * (x[1] - 4);
		dfdx[k++] = 2 * (x[2] - 1); dfdx[k++] = 2 * (x[2] - 3);
	    }
		*/
};

void Print(double *x, int n, const std::string &name = "x") {
	std::cout << name << ":";
	for (int i=0;i<n;i++) {
		std::cout << " " << x[i];
	}
	std::cout << std::endl;
}

int main(int argc, char *argv[]) {
	
	
	
	std::cout << "///////////////////////////////////////////////////" << std::endl;
	std::cout << "// Test the GCMMA algorithm" << std::endl;
	std::cout << "///////////////////////////////////////////////////" << std::endl;

	Problem toy;
	double movlim = 0.2;

	double f, fnew;
	std::vector<double> df(toy.n);
	std::vector<double> g(toy.m), gnew(toy.m);
	std::vector<double> dg(toy.n * toy.m);

	std::vector<double> x = toy.x0;
	std::vector<double> xold = x;
	std::vector<double> xnew(toy.n);


	// Print initial values
	toy.Obj(x.data(), &f, g.data());
	std::cout << "f: " << f << std::endl;
	Print(g.data(), toy.m, "g");

	// Initialize GCMMA
	GCMMASolver gcmma(toy.n, toy.m, 0, 1000, 1);
	MMASolver mma(toy.n, toy.m, 0, 1000, 1);

	double ch = 1.0;
	int maxoutit = 8;
	for (int iter = 0; ch > 0.0002 && iter < maxoutit; ++iter) {
		toy.ObjSens(x.data(), &f, g.data(), df.data(), dg.data());

		// Call the update method
		if (0) {
			// MMA version
			mma.Update(x.data(), df.data(), g.data(), dg.data(),
				toy.xmin.data(), toy.xmax.data());
		} else {
			// GCMMA version
			gcmma.OuterUpdate(xnew.data(), x.data(), f, df.data(),
				g.data(), dg.data(), toy.xmin.data(), toy.xmax.data());

			// Check conservativity
			toy.Obj(xnew.data(), &fnew, gnew.data());
			bool conserv = gcmma.ConCheck(fnew, gnew.data());
			//std::cout << conserv << std::endl;
			for (int inneriter = 0; !conserv && inneriter < 15; ++inneriter) {
				// Inner iteration update
				gcmma.InnerUpdate(xnew.data(), fnew, gnew.data(), x.data(), f,
					df.data(), g.data(), dg.data(), toy.xmin.data(), toy.xmax.data());

				// Check conservativity
				toy.Obj(xnew.data(), &fnew, gnew.data());
				conserv = gcmma.ConCheck(fnew, gnew.data());
				//std::cout << conserv << std::endl;
			}
			x = xnew;
		}

		// Compute infnorm on design change
		ch = 0.0;
		for (int i=0; i < toy.n; ++i) {
			ch = std::max(ch, std::abs(x[i] - xold[i]));
			xold[i] = x[i];
		}

		// Print to screen
		printf("it.: %d, obj.: %f, ch.: %f \n", iter, f, ch);
		Print(x.data(), toy.n);
		toy.Obj(x.data(), &f, g.data());
		std::cout << "f: " << f << std::endl;
		Print(g.data(), toy.m, "g");
		std::cout << std::endl;
	}

	return 0;
}
