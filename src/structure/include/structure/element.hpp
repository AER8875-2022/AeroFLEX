#ifndef STRUCTURE_ELEMENT_HPP
#define STRUCTURE_ELEMENT_HPP

#include "Eigen/Dense"
#include "Eigen/Geometry"
#include <vector>
#include <cmath>
#include <iostream>
#include "structure/proprietes_sections.hpp"

namespace structure {

class CBAR{

private:

public:
    PBAR ID_Prop;                    //Identifiant de prop de section
         
    Eigen::Vector3d N1;              //Coordonnées neoud #1
    Eigen::Vector3d N2;              //Coordonnées neoud #2
    Eigen::Vector3d Nmid;

    Eigen::Vector3d r1;
    Eigen::Vector3d r2;

    
    int N1_ID;                       //ID code du noeud 1
    int N2_ID;                       //ID code du noeud 2

    double L;                        //Longueur d'élément

    Eigen::VectorXd u_1;             //Déplacement du noeud 1 dans le repère local 
    Eigen::VectorXd u_2;             //Déplacement du noeud 2 dans le repère local
    
    Eigen::Vector3d V;               //Vecteur d'orientation
    Eigen::VectorXd u_mid;           //Déplacement du noeud moyen
    
    Eigen::Quaterniond q_1;          //Quaternion du noeud 1 
    Eigen::Quaterniond q_1_rot_prime;

    Eigen::Quaterniond q_mid;

    Eigen::Quaterniond q_2;          //Quaternion du noeud 2 
    Eigen::Quaterniond q_2_rot_prime;

    Eigen::Quaterniond q_e;          //Quaternion initiale de l'élément

    Eigen::MatrixXd K_elem_local;    //Matrice de rigidité dans le repère local
    Eigen::MatrixXd K_elem_global;   //Matrice de rigidité dans le repère global
    Eigen::MatrixXd T_Rotation;      //Matrice de rotation pour la matrice de rigidité

    CBAR(){}

    CBAR(PBAR id_prop, Eigen::Vector3d noeud1, Eigen::Vector3d noeud2,Eigen::VectorXd v, int n1_id, int n2_id){
        ID_Prop = id_prop;
        N1      = noeud1;
        N2      = noeud2;
        Nmid    = 0.5*(N1+N2);
        r1      = N1 - Nmid;
        r2      = N2 - Nmid;
        V       = v;
        N1_ID   = n1_id;
        N2_ID   = n2_id;
        L       = sqrtl( (noeud1(0)-noeud2(0))*(noeud1(0)-noeud2(0))+   (noeud1(1)-noeud2(1))*(noeud1(1)-noeud2(1)) + (noeud1(2)-noeud2(2))*(noeud1(2)-noeud2(2)) );
        
        set_K_local();
        set_T_Rotation();
        set_K_elem_global();
        
        u_1 = Eigen::VectorXd::Zero(6);
        u_2 = Eigen::VectorXd::Zero(6);
        
    };


    void set_K_local()
    {
        const double L_inv= 1/L;

        const double X = ID_Prop.A * ID_Prop.E * L_inv;

        const double Y2 = 6. * ID_Prop.E * ID_Prop.Iz * L_inv * L_inv;
        const double Y1 = Y2 * 2. * L_inv;
        const double Y3 = 4.  * ID_Prop.E * ID_Prop.Iz * L_inv;
        const double Y4 = Y3/2. ;
    
        const double Z2 = 6. * ID_Prop.E * ID_Prop.Iy * L_inv * L_inv;
        const double Z1 = Z2 * 2. * L_inv;
        const double Z3 = 4. * ID_Prop.E * ID_Prop.Iy * L_inv;
        const double Z4 = Z3/2. ;

        const double S = ID_Prop.G * ID_Prop.J * L_inv;

        K_elem_local=Eigen::MatrixXd::Zero(12,12);
        
        K_elem_local << X,  0,  0, 0,  0,  0,-X,  0,  0,  0,  0,  0,
                        0, Y1,  0, 0,  0, Y2, 0,-Y1,  0,  0,  0, Y2,
                        0,  0, Z1, 0,-Z2,  0, 0,  0,-Z1,  0,-Z2,  0,
                        0,  0,  0, S,  0,  0, 0,  0,  0, -S,  0,  0,
                        0,  0,-Z2, 0, Z3,  0, 0,  0, Z2,  0, Z4,  0,
                        0, Y2,  0, 0,  0, Y3, 0,-Y2,  0,  0,  0, Y4,
                       -X,  0,  0, 0,  0,  0, X,  0,  0,  0,  0,  0,
                        0,-Y1,  0, 0,  0,-Y2, 0, Y1,  0,  0,  0,-Y2,
                        0,  0,-Z1, 0, Z2,  0, 0,  0, Z1,  0, Z2,  0,
                        0,  0,  0,-S,  0,  0, 0,  0,  0,  S,  0,  0,
                        0,  0,-Z2, 0, Z4,  0, 0,  0, Z2,  0, Z3,  0,
                        0, Y2,  0, 0,  0, Y4, 0,-Y2,  0,  0,  0, Y3;
    };
    
    void set_T_Rotation()
    {
        const double L_inv= 1/L;
        Eigen::MatrixXd Lambda= Eigen::MatrixXd::Zero(3,3);

        //Vecteur X_bar , soit le X local (Unitaire)
        const Eigen::Vector3d unit_N1_to_N2 = L_inv*(N2 - N1);         
        const double l1 = unit_N1_to_N2(0);
        const double m1 = unit_N1_to_N2(1);
        const double n1 = unit_N1_to_N2(2);    

        //Vecteur Y_bar est parallèle au vecteur V. Juste à le rendre unitaire
        V= V/V.norm();
        const double l2 = V(0);
        const double m2 = V(1);
        const double n2 = V(2);

        //Vecteur Z_bar est ortogonal au X_bar et Y_bar. Pour trouver ses composantes, faire un produit vectoriel.
        const Eigen::Vector3d cross_product = unit_N1_to_N2.cross(V);
        const double l3 = cross_product(0);
        const double m3 = cross_product(1);
        const double n3 = cross_product(2);      

        Lambda   << l1, m1, n1,
                    l2, m2, n2,  
                    l3, m3, n3;


        T_Rotation                   = Eigen::MatrixXd::Zero(12,12);
        T_Rotation.block(0, 0, 3, 3) = Lambda;
        T_Rotation.block(3, 3, 3, 3) = Lambda;
        T_Rotation.block(6, 6, 3, 3) = Lambda;
        T_Rotation.block(9, 9, 3, 3) = Lambda;

        set_qe_From_Rotation_Matrix(Lambda);            
    };


    void set_K_elem_global()  //Définir la matrice de rigidité de l'élément dans le repère global
    {   
        K_elem_global = T_Rotation.transpose() * K_elem_local * T_Rotation ;
        
    };



    void set_u_i(Eigen::VectorXd Delta_dep1 , Eigen::VectorXd Delta_dep2 )  //
    {
        u_1   += Delta_dep1;  //Déplacement du noeud 1 dans le repère global
        u_2   += Delta_dep2;  //Déplacement du noeud 2 dans le repère global
        u_mid  = 0.5 * (u_1 + u_2);  //Déplacement du noeud milieu dans le repère global
    };


    void set_qe_From_Rotation_Matrix(Eigen::Matrix3d mat)  //Set q_e
    {
        const double trace = mat(0,0) + mat(1,1) + mat(2,2);
        double cst, s, v_x, v_y, v_z;

        if (trace > 0.0)
        {
            cst = 0.5  / sqrt(trace + 1.0);
            s   = 0.25 / cst;

            v_x = (mat(1,2) - mat(2,1)) * cst;
            v_y = (mat(2,0) - mat(0,2)) * cst;
            v_z = (mat(0,1) - mat(1,0)) * cst;
        
        }
        else
        {
            if (mat(0,0) > mat(1,1) && mat(0,0) > mat(2,2))
            {
                cst = 2.0 * sqrt(1.0 + mat(0,0) - mat(1,1) - mat(2,2));
                s   = ( mat(2,1) - mat(1,2) ) / cst;

                v_x = 0.25 * cst ;
                v_y = ( mat(0,1) + mat(1,0) ) / cst;
                v_z = ( mat(0,2) + mat(2,0) ) / cst;
            }
            else if(mat(1,1) > mat(2,2))
            {
                cst = 2.0 * sqrt( 1.0 + mat(1,1) - mat(0,0) - mat(2,2) );
                s   = ( mat(0,2) - mat(2,0) ) / cst;
                v_x = ( mat(0,1) + mat(1,0) ) / cst;
                v_y = 0.25 * cst;
                v_z = ( mat(1,2) + mat(2,1) ) / cst;
            }
            else
            {
                cst = 2.0 * sqrt(1.0 + mat(2,2) - mat(0,0)- mat(1,1));
                s   = ( mat(1,0) - mat(0,1) ) / cst;
                v_x = ( mat(0,2) + mat(2,0) ) / cst;
                v_y = ( mat(1,2) + mat(2,1) ) / cst;
                v_z = 0.25 * cst;
            }
        }
        
        q_e = Eigen::Quaterniond(s, v_x, v_y, v_z);
        q_1 = q_e;
        
        q_2 = q_e;
        q_mid = q_e;

        
        
    }

    void set_q1_And_q2(Eigen::Quaterniond delta_q_1 ,Eigen::Quaterniond delta_q_2 )
    {
        q_1 = delta_q_1 * q_1;
        q_2 = delta_q_2 * q_2;
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
    

   void set_qmid_From_Interpolation()  //Trouver q_mid à partir de q_1 et q_2
    {
        if ( abs(q_1.w()-q_2.w()) <1e-14 && abs(q_1.x()-q_2.x()) < 1e-14 && abs(q_1.y()-q_2.y()) < 1e-14 && abs(q_1.z()-q_2.z()) < 1e-14)
        {
            q_mid = q_1;
        }
        else{
            
            double theta = acos((q_1.w()*q_2.w() + q_1.x()*q_2.x() + q_1.y()*q_2.y() + q_1.z()*q_2.z())/(q_1.norm()*q_2.norm()));
                    
            if (abs(theta) > 1e-14)
            {   
                Eigen::Quaterniond p;
                const double cst   = sin(0.5*theta)/sin(theta);
                
                p.w() = cst*q_1.w() + cst*q_2.w();
                p.x() = cst*q_1.x() + cst*q_2.x();
                p.y() = cst*q_1.y() + cst*q_2.y();
                p.z() = cst*q_1.z() + cst*q_2.z();

                q_mid =p;
            }
            else
            {   
                q_mid =  q_1;
            } 
        }
    };

    void set_Quaternion_Local_Rotations() //Puisqu'on utilise des quaternion unitaire, le conjugate et la transpose sont équivalent.
    {   
        q_1_rot_prime = q_mid.conjugate()*q_e.conjugate() * q_1 * q_e;
        q_2_rot_prime = q_mid.conjugate()*q_e.conjugate() * q_2 * q_e;
        
    };

    Eigen::Vector3d get_Euler_Angles_From_Local_Rotation(Eigen::Quaterniond q_n_rot_local)  //Trouve les déformations angulaires à partir des q_rot,i'
    {
        const double s = q_n_rot_local.w();
        const double a = q_n_rot_local.x();
        const double b = q_n_rot_local.y();
        const double c = q_n_rot_local.z();

        Eigen::Vector3d THETA_local;

        THETA_local(0) = atan2(2 * s * a + 2 * b * c , 1 - 2 * ( a*a + b*b) );               
        THETA_local(1) = asin( 2 * s * b - 2 * c * a);
        THETA_local(2) = atan2(2 * s * c + 2 * a * b , 1 - 2 * ( b*b + c*c));

        return THETA_local;
    };

    Eigen::VectorXd get_Deformation_Local_Ref()
    {    
        const Eigen::Vector3d theta_1_local = get_Euler_Angles_From_Local_Rotation(q_1_rot_prime);
        const Eigen::Vector3d theta_2_local = get_Euler_Angles_From_Local_Rotation(q_2_rot_prime);

        

        const Eigen::Matrix3d R_ec          = get_Rotation_Matrix_From_Quaternion(q_mid * q_e.inverse());              //R(q_mid)
        
        const Eigen::Vector3d dr_1          = R_ec*r1 - r1;                                                           //Déplacement induit par la rotation
        const Eigen::Vector3d dr_2          = R_ec*r2 - r2;                                                           //Déplacement induit par la rotation 
    
        const Eigen::Matrix3d R_gc          = get_Rotation_Matrix_From_Quaternion(q_mid);

        

        const Eigen::Vector3d d_1_prime     =  R_gc.transpose() * (u_1.segment(0,3) - u_mid.segment(0,3) - dr_1);     //Déplacement du noeud causé pas les déformations dans le repère de l'élément
        const Eigen::Vector3d d_2_prime     =  R_gc.transpose() * (u_2.segment(0,3) - u_mid.segment(0,3) - dr_2);     //Déplacement du noeud causé pas les déformations dans le repère de l'élément


        // std::cout<<"============\n"<<std::endl;
        // std::cout<<dr_1.transpose()<<std::endl;
        // std::cout<<dr_2.transpose()<<std::endl;
        // std::cout<<d_1_prime.transpose()<<std::endl;
        // std::cout<<d_2_prime.transpose()<<std::endl;

        Eigen::VectorXd D_local(12);
        
        D_local.segment(0,3) = d_1_prime;
        D_local.segment(3,3) = theta_1_local;
        D_local.segment(6,3) = d_2_prime;
        D_local.segment(9,3) = theta_2_local;

        return D_local;
    };

    Eigen::VectorXd get_Force_In_GlobalRef(Eigen::VectorXd D_local)
    {   
        
        Eigen::VectorXd F_int_elem_local = K_elem_local * D_local;
                    
        Eigen::MatrixXd Rot  = Eigen::MatrixXd::Zero(12,12);
        Eigen::Matrix3d Diag = get_Rotation_Matrix_From_Quaternion(q_mid);
        

        Rot.block(0, 0, 3, 3) = Diag;
        Rot.block(3, 3, 3, 3) = Diag;
        Rot.block(6, 6, 3, 3) = Diag;
        Rot.block(9, 9, 3, 3) = Diag;

        return Rot*F_int_elem_local;  
    }
        
};

}

#endif
