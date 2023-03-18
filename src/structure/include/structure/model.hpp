#ifndef STRUCTURE_MODELE_HPP
#define STRUCTURE_MODELE_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <vector>
#include <cmath>
#include <map>
#include <cstdlib>
#include <structure/element.hpp>
#include <structure/spc.hpp>
#include <structure/loads.hpp>
#include <fstream>
#include <Eigen/SparseLU>
#include <atomic>
#include <thread>
#include "common_aeroflex.hpp"

namespace structure{

class MODEL{

public:
    std::map<int,int>                     indexation_switch;      //  user_id(ni) --> code_id(numéro du ddl)
    std::map<int, Eigen::Vector3d>                 Grid_MAP;      //  user_id(ni) --> Initial position(x,y,z)
    std::map<int, MAT1>                            MAT1_MAP;      //  MAT1_ID --> MAT1
    std::map<int, PBAR>                            PBAR_MAP;      //  PBAR_ID --> PBAR
    std::map<int, CBAR>                            CBAR_MAP;      //  CBAR_ID --> CBAR
    std::map<int, Eigen::Quaterniond>        QUATERNION_MAP;      //  code_id --> CBAR 
    std::vector<SPC1>                             SPC1_LIST;      //   
    std::vector<FORCE>                           FORCE_LIST;      //
    std::vector<MOMENT>                         MOMENT_LIST;      //

    int Nbr_Element;                        //Nombre d'éléments
    unsigned int Nbr_Noeud;                 //Nombre de noeud 

    Eigen::VectorXd                          Forces;
    Eigen::SparseMatrix<double>     K_Global_sparse;
    Eigen::SparseMatrix<double>      K_Final_sparse;
    Eigen::VectorXd                    Displacement;

    // GUI Handlers
    GUIHandler &gui;
    std::atomic<int> &iters;
    std::vector<double> &residuals;

    MODEL(GUIHandler &gui, std::atomic<int> &iters, std::vector<double> &residuals)
        : gui(gui), iters(iters), residuals(residuals) {};

    MODEL(std::string namefile, GUIHandler &gui, std::atomic<int> &iters, std::vector<double> &residuals)
        : gui(gui), iters(iters), residuals(residuals) {
        
        read_data_file(namefile);
        set_K_global();
        set_Load_Vector_From_Load_Objects();    
        set_K_Final_sparse();
    };

    struct PBAR_param{
        int PBAR_id;
        int MAT1_id;
        
        Eigen::VectorXd PARAM;  // a, iz, iy, j,
    };

    struct CBAR_param{
        int             CBAR_id;
        int             PBAR_id;
        int             N1_user_ID;
        int             N2_user_ID;
        Eigen::Vector3d V;

    };

    struct SPC1_param{
        int    User_node_id;   //Utilisateur
        std::string    code;   //
    };

    struct FORCE_param{
        int             user_node_id;         
        double          Norm;              
        Eigen::Vector3d Direction;
    };

    struct MOMENT_param{
        int             user_node_id;         
        double          Norm;              
        Eigen::Vector3d Direction;
    };

    void set_Load_Vector_From_Vector(Eigen::VectorXd New_F)
    {
        if(Forces.size() == New_F.size())
        {
            Forces = New_F;
        }
    }

    
    void set_K_global(){

        K_Global_sparse = Eigen::SparseMatrix<double>( Nbr_Noeud * 6, Nbr_Noeud * 6 );
        for( auto& [cbar_id, elem] : CBAR_MAP)
        {
            int n1  = elem.N1_ID; 
            int n2  = elem.N2_ID;
            Eigen::Matrix3d Diag  = elem.get_Rotation_Matrix_From_Quaternion(elem.q_mid);
            Eigen::MatrixXd Rot=Eigen::MatrixXd::Zero(12,12);

            Rot.block(0,0,3,3) = Diag;
            Rot.block(3,3,3,3) = Diag;
            Rot.block(6,6,3,3) = Diag;
            Rot.block(9,9,3,3) = Diag;

            Eigen::MatrixXd Tempo = Rot * elem.K_elem_local * Rot.transpose();
            
            for (unsigned int j = 0; j < 6; j++)
            {  
                for (unsigned int k = 0 ; k < 6 ; k++)
                {
                    K_Global_sparse.coeffRef(6*n1+j,6*n1+k) += Tempo(j,k);
                    K_Global_sparse.coeffRef(6*n2+j,6*n2+k) += Tempo(6+j,6+k);
                    K_Global_sparse.coeffRef(6*n1+j,6*n2+k) += Tempo(j,6+k);
                    K_Global_sparse.coeffRef(6*n2+j,6*n1+k) += Tempo(6+j,k);
                }
            }                 
        }
    }

    void set_Load_Vector_From_Load_Objects()
    {      
        Forces.setZero(6*Nbr_Noeud);
        // #pragma omp parallel for
        for(int i=0 ; i < FORCE_LIST.size();i++)
        { 
            FORCE f_obj       = FORCE_LIST[i];
            Eigen::Vector3d f = f_obj.get_xyz_force();
            Forces.segment(f_obj.Node_ID *6,3) += f;
        };
        // #pragma omp parallel for
        for(int i=0 ; i<MOMENT_LIST.size();i++)
        {
            MOMENT m_obj       = MOMENT_LIST[i];
            Eigen::Vector3d m  = m_obj.get_xyz_moment();
            Forces.segment((m_obj.Node_ID)*6+3 ,3) += m;
        };
    };

    Eigen::VectorXd apply_SPC1_Forces(Eigen::VectorXd force)
    {   
        // #pragma omp parallel for
        for (int i=0 ; i<SPC1_LIST.size();i++)
        {   
            SPC1 spc_obj     = SPC1_LIST[i];
            std::string code = spc_obj.CODE;
            for(char j : code)
            {   
                int J = j - '0';
                int ddl      = (6*spc_obj.Node_ID) + J - 1;
                force(ddl)   = 0;
            }
        }
        return force;
    }

    void set_K_Final_sparse()
    {
        K_Final_sparse = K_Global_sparse;
        // #pragma omp parallel for
        for (int i=0 ; i<SPC1_LIST.size();i++)
        {
            SPC1 spc_obj     = SPC1_LIST[i];
            
            std::string code = spc_obj.CODE;
            for(char j : code)
            {   
                int J = j - '0';
                int ddl      = (6*spc_obj.Node_ID) +J -1;
                Forces(ddl)  = 0;
                for (int k = 0; k < K_Final_sparse.cols(); ++k)
                {
                    K_Final_sparse.coeffRef(ddl,k) = 0;
                    K_Final_sparse.coeffRef(k,ddl) = 0;
                }      
                K_Final_sparse.coeffRef(ddl,ddl)   = 1;
            }
        }
        K_Final_sparse.prune(0.0);
    }

    void Rotate_Ext_Loads()
    {   
        for (unsigned int i = 0; i < Nbr_Noeud; i++)
        {
            Eigen::Matrix3d diag  = get_Rotation_Matrix_From_Quaternion(QUATERNION_MAP[i]); 
            
            Eigen::MatrixXd Rot   = Eigen::MatrixXd::Zero(6,6);
            Rot.block(0,0,3,3)    = diag;
            Rot.block(3,3,3,3)    = diag;
            
            Forces.segment(i*6,6) = Rot * Forces.segment(i*6,6);            
        }
    }

    Eigen::Matrix3d get_Rotation_Matrix_From_Quaternion(Eigen::Quaterniond q)
    {
        const double s = q.w();
        const double x = q.x();
        const double y = q.y();
        const double z = q.z();

        Eigen::Matrix3d Matrix = Eigen::Matrix3d::Zero(3,3);

        Matrix << 1. - 2.*(y*y + z*z), 2.*(x*y - s*z)     , 2.*(x*z + s*y),
                  2.*(x*y + s*z)     , 1. - 2.*(x*x + z*z), 2.*(y*z - s*x),
                  2.*(x*z - s*y)     , 2.*(y*z + s*x)     , 1. - 2.*(x*x + y*y);

        return Matrix;
    };


    Eigen::VectorXd get_Solve(Eigen::VectorXd f)
    {   
        Eigen::SparseLU<Eigen::SparseMatrix<double>> Solver;
        Eigen::SparseMatrix<double> Tempo = K_Final_sparse; 
        Solver.compute(Tempo);
        
        return Solver.solve(f); 
    }

    Eigen::VectorXd get_Lin_Solve()
    {   
        return get_Solve(Forces); 
    }
    
    Eigen::VectorXd get_NonLin_Solve(int Max_load_step, double tol, double amor)
    {
        Eigen::VectorXd Dep = Eigen::VectorXd(6 * Nbr_Noeud);
        Eigen::VectorXd Forces_int(6 * Nbr_Noeud);
        Eigen::VectorXd Forces_diff(6 * Nbr_Noeud );

        //Dep.setZero();
        Forces_int.setZero();
        Forces_int.setZero();

        std::vector<int> CBAR_keys;
        CBAR_keys.reserve(CBAR_MAP.size());
        for (auto& element : CBAR_MAP) 
        {
                CBAR_keys.push_back(element.first);
        }

        for (int i = 0; i < Nbr_Noeud; i++)
        {
            Eigen::Quaterniond delta_q(0.0,0.0,0.0,0.0);
            QUATERNION_MAP[i] = delta_q ;
        }
            

        for (double Load_Step = 1.; Load_Step <= Max_load_step; Load_Step ++)
        {
            std::cout<<"\n============"<<std::endl;
            std::cout<<"Load step: "<< Load_Step <<" / "<< Max_load_step <<std::endl;
            Forces_diff = (Load_Step/Max_load_step)*(Forces) - Forces_int;
            
            set_K_global();
            set_K_Final_sparse();
            Eigen::VectorXd Delta_dep_full; 
            Eigen::VectorXd Delta_dep_amor = get_Solve(Forces_diff); 
              
            double Residu = 1.0;    

            // Main solving loop
            do {
                while (gui.signal.pause) std::this_thread::sleep_for(std::chrono::milliseconds(100));
                
                //Delta déplacements
                Dep += Delta_dep_amor; 
                std::cout<<"==============="<<std::endl;
                std::cout<<Dep.tail(3).transpose()<<std::endl;
                set_Quaternion_Map(Dep); 

                //#pragma omp parallel for
                for (int i = 0; i < CBAR_keys.size(); ++i)
                {   
                    int key     = CBAR_keys[i];
                    CBAR& value = CBAR_MAP[key]; 

                    int n1,n2;
                    n1 = value.N1_ID;
                    n2 = value.N2_ID;
                    
                    Eigen::VectorXd delta_dep1 = Delta_dep_amor.segment(n1*6,6);
                    Eigen::VectorXd delta_dep2 = Delta_dep_amor.segment(n2*6,6);  
                     
                                   
                    value.set_u_i(delta_dep1,delta_dep2);                                      //Set u_1 and u_2
                    value.set_q1_And_q2(QUATERNION_MAP[n1],QUATERNION_MAP[n2]);                //Set q_1 and q_2                    
                    value.set_qmid_From_Interpolation();                                       //Set q_mid                                         
                    value.set_Quaternion_Local_Rotations();                                    //Set q_1_rot_prime and q_2_rot_prime 
                    
                    Eigen::VectorXd F_elem_global_ref = value.get_Force_In_GlobalRef();
                    
                    #pragma omp critical
                    Forces_int.segment(6*n1,6)       += F_elem_global_ref.segment(0,6);      
                    #pragma omp critical
                    Forces_int.segment(6*n2,6)       += F_elem_global_ref.segment(6,6);   
                }
                
                Forces_int = apply_SPC1_Forces(Forces_int);
                Rotate_Ext_Loads();
                Forces_diff = (Load_Step/Max_load_step)*(Forces) - Forces_int;
                
                set_K_global();
                set_K_Final_sparse();
                Delta_dep_full = get_Solve(Forces_diff); 
                Residu = std::sqrt(Delta_dep_full.transpose()*Delta_dep_full);
                Forces_int.setZero();
                
                if (Residu < 100.*tol) Delta_dep_amor = amor*Delta_dep_full;
                else Delta_dep_amor = Delta_dep_full;

                if (iters%100 == 0 || Residu < tol){
                std::cout << "Iteration " << iters << std::endl;
                std::cout << "\t Residual = " << Residu << std::endl;
                };
                
                // Incrementing current iteration
                iters++;
                residuals.push_back(Residu);
            } while((Residu> tol) && (!gui.signal.stop));
        }
        return Dep;
    }

    void set_Quaternion_Map(Eigen::VectorXd delta_dep){

        //#pragma omp parallel for
        for (int i = 0; i < Nbr_Noeud; i++)
        {   
            double rx = delta_dep(i*6 +3);
            double ry = delta_dep(i*6 +4);
            double rz = delta_dep(i*6 +5);
            
            double vx = cos(0.5*rz) * cos(0.5*ry) * sin(0.5*rx) - sin(0.5*rz) * sin(0.5*ry) * cos(0.5*rx);
            //if(abs(vx)<1E-12) vx=0.0;

            double vy = cos(0.5*rz) * sin(0.5*ry) * cos(0.5*rx) + sin(0.5*rz) * cos(0.5*ry) * sin(0.5*rx);
            //if(abs(vy)<1E-12) vy=0.0;

            double vz = sin(0.5*rz) * cos(0.5*ry) * cos(0.5*rx) - cos(0.5*rz) * sin(0.5*ry) * sin(0.5*rx);
            //if(abs(vz)<1E-12) vz=0.0;

            double s  = cos(0.5*rz) * cos(0.5*ry) * cos(0.5*rx) + sin(0.5*rz) * sin(0.5*ry) * sin(0.5*rx);
            //if(abs(s)<1E-12) s=0.0;

            Eigen::Quaterniond delta_q(s,vx,vy,vz);

            if(i== (Nbr_Noeud-1))
            {   
                std::cout<<delta_q<<std::endl;
            }
            QUATERNION_MAP[i] = delta_q ;
        }
    }
    
    void read_data_file(std::string namefile)
    {   
        std::ifstream file(namefile);
        if (!file.is_open()) 
        {
            std::cerr << "\n\033[1;31m ->FEM ERROR: name file \"" << namefile << "\" not found! \033[0m" << std::endl;
            exit(1);
        }
        
        //Delimiter
        char delimiter = ',';
        
        //Seperate variables
        std::vector<std::string> line;

        //Initialisation
        Nbr_Noeud    =0;
        Nbr_Element  =0;
        
        //Initialisation des vecteur pour stocker les paramètres et les créer à la fin.
        std::vector<PBAR_param>   PBAR_stock;
        std::vector<CBAR_param>   CBAR_stock;
        std::vector<SPC1_param>   SPC1_stock;
        std::vector<FORCE_param>  FORCE_stock;
        std::vector<MOMENT_param> MOMENT_stock;

        //Line
        std::string fileline;
        while (std::getline(file, fileline)) 
        {
            if (!fileline.size())
            {
                continue;
            }

            //Clear the vector
            line.clear(); 

            //Delete de spaces
            fileline.erase(std::remove(fileline.begin(),fileline.end(),' '),fileline.end());
            
            std::string word;

            std::stringstream lineStream(fileline);
            while (std::getline(lineStream, word, delimiter)) 
            {
                line.push_back(word);
            }
            
            if( line[0]=="GRID")
            {   
                int data_id = std::stoi(line[1]);
                double             x = std::stod(line[3]);
                double             y = std::stod(line[4]);
                double             z = std::stod(line[5]);            
                
                Eigen::Vector3d node(x,y,z);

                //Create a node
                Grid_MAP[data_id] = node;
                indexation_switch[data_id] = Nbr_Noeud;
                Nbr_Noeud+=1;
            }
            if(line[0]=="MAT1")
            {
                int data_id = std::stoi(line[1]);
                double             e = std::stod(line[2]);
                double             g = std::stod(line[3]);
                
                 
                //Create a MAT1
                MAT1_MAP[data_id] = MAT1(e,g);
            }
            if(line[0]=="PBAR")
            {   
                int data_id = std::stoi(line[1]);
                int MAT1_id = std::stoi(line[2]);
                

                double a             = std::stod(line[3]);
                double iz            = std::stod(line[4]);
                double iy            = std::stod(line[5]);
                double j             = std::stod(line[6]);
                Eigen::VectorXd param(4);

                param(0) = a;
                param(1) = iz;
                param(2) = iy;
                param(3) = j;

                PBAR_param p = {data_id,MAT1_id,param};

                PBAR_stock.push_back(p);
            }

            if(line[0]=="CBAR")
            {   
                 
                int data_id    = std::stoi(line[1]);
                
                int PBAR_id    = std::stoi(line[2]);
                int n1_user_id = std::stoi(line[3]);
                int n2_user_id = std::stoi(line[4]);
                
                double x       = std::stod(line[5]);
                double y       = std::stod(line[6]);
                double z       = std::stod(line[7]);
                Eigen::Vector3d  V(x,y,z);

                CBAR_param c = {data_id, PBAR_id, n1_user_id,n2_user_id,V};

                CBAR_stock.push_back(c);
                
            }
            if(line[0]=="SPC1")
            {   
                int data_id           = std::stoi(line[1]);
                std::string CODE      = line[2];
                
                for ( int i = 3; i < line.size(); i++)
                {   
                    
                    if (!(static_cast<int>(line[i][0])==13))
                    {
                        int node_id = std::stoi(line[i]);   //Node utilisateur imposé à zéro
                        SPC1_param spc1 = { node_id , CODE};
                        SPC1_stock.push_back(spc1);
                    }
                }             
            }
            if(line[0]=="FORCE")
            {   
                int data_id               = std::stoi(line[1]);
                int user_node_id          = std::stoi(line[2]);
                double norm               = std::stod(line[4]);

                double x                  = std::stod(line[5]);
                double y                  = std::stod(line[6]);
                double z                  = std::stod(line[7]);
                Eigen::Vector3d direction(x,y,z);

                FORCE_param f = {user_node_id,norm,direction};
                FORCE_stock.push_back(f);           
            }
            if(line[0]=="MOMENT")
            {   
                int data_id               = std::stoi(line[1]);
                int user_node_id          = std::stoi(line[2]);
                double norm               = std::stod(line[4]);

                double x                  = std::stod(line[5]);
                double y                  = std::stod(line[6]);
                double z                  = std::stod(line[7]);
                Eigen::Vector3d direction(x,y,z);

                MOMENT_param m = {user_node_id,norm,direction}; 
                MOMENT_stock.push_back(m);              
            }

        }
        file.close();
        //Once the file is over:

        //Create PBAR
        
        for (int i = 0; i < PBAR_stock.size(); i++)
        {   
            int pbar_id = PBAR_stock[i].PBAR_id;
            int mat1_id = PBAR_stock[i].MAT1_id;
            Eigen::VectorXd param = PBAR_stock[i].PARAM;  // a, iz, iy, j,
            
            //Create the objet
            PBAR Pi = PBAR(param(0), param(1), param(2), param(3),MAT1_MAP[mat1_id]);
            
            //Store the objet in the map
            PBAR_MAP[pbar_id] = Pi;
        }
        
        //Create CBAR
        for (int i = 0; i < CBAR_stock.size(); i++)
        {
            int cbar_id        = CBAR_stock[i].CBAR_id;
            int pbar_id        = CBAR_stock[i].PBAR_id;
            int n1_user_id     = CBAR_stock[i].N1_user_ID;
            int n2_user_id     = CBAR_stock[i].N2_user_ID;
            Eigen::Vector3d v  = CBAR_stock[i].V;

            
            int n1_code_id     = indexation_switch[n1_user_id];
            int n2_code_id     = indexation_switch[n2_user_id];

            Eigen::Vector3d n1 = Grid_MAP[n1_user_id];
            Eigen::Vector3d n2 = Grid_MAP[n2_user_id]; 

            PBAR p_use         = PBAR_MAP[pbar_id];

            //Create the objet
            CBAR Ci = CBAR(p_use,n1,n2,v,n1_code_id,n2_code_id);

            //Store the objet in the map
            CBAR_MAP[cbar_id] = Ci;
            Nbr_Element++;
        }

        //Create SPC1
        for (int i = 0; i < SPC1_stock.size(); i++)
        {
            int user_node_id = SPC1_stock[i].User_node_id;
            std::string code = SPC1_stock[i].code;

            //Get the position of the node in the matrix/force vector
            int code_node_id = indexation_switch[user_node_id];
            
            //Create the objet
            SPC1 spc1 = SPC1(code_node_id,code);

            SPC1_LIST.push_back(spc1);
        }
        //Create FORCE
        for (int i = 0; i < FORCE_stock.size(); i++)
        {
            int user_node_id  = FORCE_stock[i].user_node_id;
            double norm       = FORCE_stock[i].Norm;
            Eigen::Vector3d v = FORCE_stock[i].Direction;

            //Get the position of the node in the matrix/force vector
            int code_node_id = indexation_switch[user_node_id];
            
            //Create the objet
            FORCE f = FORCE(code_node_id,norm,v);

            FORCE_LIST.push_back(f);
        }
        //Create MOMENT
        for (int i = 0; i < MOMENT_stock.size(); i++)
        {
            int user_node_id  = MOMENT_stock[i].user_node_id;
            double norm       = MOMENT_stock[i].Norm;
            Eigen::Vector3d v = MOMENT_stock[i].Direction;

            //Get the position of the node in the matrix/force vector
            int code_node_id = indexation_switch[user_node_id];
            
            //Create the objet
            MOMENT m = MOMENT(code_node_id,norm,v);

            MOMENT_LIST.push_back(m);
        }
        
    }
};

}

#endif
