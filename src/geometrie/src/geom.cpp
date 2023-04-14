#include "geometrie/geom.hpp"
#include <cmath>
#include <iostream>
#include <ostream>
#include <string.h>
#include <vector>
#include <sstream>


#ifndef M_PI
#define M_PI 3.141592653589793115997963468544185161590576171875
#endif

using namespace geom;

//std::string body_type = "General";
Geom::Geom(GUIHandler &gui): gui(gui), WING_RIGHT("General"), WING_LEFT("General") {}

void Geom::Geom_gen() {
    // convertion angles - distances
    double xle;
    if (settings.fleche==0.0){
        xle = 0.0;
    } else {
        xle = settings.envergure/(2*tan(settings.fleche*(M_PI/180)));
    }
    double zn =settings.envergure*(sin(settings.dihedre*(M_PI/180))/2);
    double d_alpha = settings.twist*(M_PI/180);

    //Define Geometry
    // std::string body_type = "General";
    // Body WING_RIGHT(body_type);
    // Body WING_LEFT(body_type);

    int nb_profils = 1; //to be changed for future work
    //double perc_span = 1.0;
    for (int k=0; k<nb_profils; k++){
        if(settings.S_type == 1){
            //WR
            //#1
            WING_RIGHT.add_wing_surface_cst(
            100,                            
            20,
            settings.envergure,             
            0,                              
            {settings.cr,settings.ct},      
            {settings.z_te},                         
            {settings.r_le},                        
            {settings.Beta},                         
            {0.0,xle},                                           
            {0.0,zn},                          
            {0.0,d_alpha}                       
            );
            // //#2
            // WING_RIGHT.add_wing_surface_cst(
            // 100,                            
            // 5,
            // settings.envergure,             
            // settings.envergure,                              
            // {settings.cr,settings.ct},      
            // {settings.z_te},                         
            // {settings.r_le},                        
            // {settings.Beta},                         
            // {0.0,xle},                                           
            // {0.0,zn},                          
            // {0.0,d_alpha}                       
            // );

            //WL
            //#1
            WING_LEFT.add_wing_surface_cst(
            100,                            
            20,
            settings.envergure,             
            0,                              
            {settings.cr,settings.ct},      
            {settings.z_te},                         
            {settings.r_le},                        
            {settings.Beta},                         
            {0.0,xle},                                           
            {0.0,zn},                          
            {0.0,d_alpha}                        
            );
            // //#2
            // WING_LEFT.add_wing_surface_cst(
            // 100,                            
            // 5,
            // settings.envergure,             
            // settings.envergure,                              
            // {settings.cr,settings.ct},      
            // {settings.z_te},                         
            // {settings.r_le},                        
            // {settings.Beta},                         
            // {0.0,xle},                                           
            // {0.0,zn},                          
            // {0.0,d_alpha}                        
            // );
            profile_name.push_back("CST");
        } else {
            //WR
            WING_RIGHT.add_wing_surface_naca(
            100,                            
            20,                              
            settings.envergure,             
            0, 
            {settings.cr,settings.ct}, 
            {settings.c},                   
            {settings.p},                        
            {settings.t},                     
            {0.0,xle},                                           
            {0.0,zn},                          
            {0.0,d_alpha}                      
            );

            //WL
            WING_LEFT.add_wing_surface_naca(
            100,                            
            20,                              
            settings.envergure,             
            0, 
            {settings.cr,settings.ct}, 
            {settings.c},                   
            {settings.p},                        
            {settings.t},                     
            {0.0,xle},                                           
            {0.0,zn},                          
            {0.0,d_alpha}                     
            );

            std::ostringstream ss1,ss2,ss3;
            ss1 << settings.c;
            ss2 << settings.p;
            ss3 << settings.t;
            ss1.precision(0);
            ss2.precision(0);
            ss3.precision(0);
            profile_name.push_back("NACA"+ss1.str()+ss2.str()+ss3.str());
            std::cout<<"NACA Name generated"<<std::endl;

        }
    }
    WING_LEFT.mirror_body(); 
    gui.msg.push("[GEOM] Geometry generated");
}


void Geom::Geom_mesh(bool viscous) {
    //VLM ------------- //Needs to be ("cartesian", "cartesian")
        // Initiate and build Mesh object
    Mesh mesh = Mesh();
    get_surface_mesh_tables(WING_RIGHT, mesh, "Structured", true);
    get_surface_mesh_tables(WING_LEFT, mesh, "Structured", true);
        // Cleanup Mesh
    double tol = 1e-14;
    int nb_subdivisions = 10;
    mesh.stitch_used_vertices(tol, nb_subdivisions);
    mesh.delete_unused_vertices();
    mesh.cleanup_doublets();
        // Export mesh object to .dat file
    export_mesh("mesh_vlm.dat", mesh);
    gui.msg.push("[GEOM] Mesh for VLM solver generated");
    gui.msg.push("[GEOM] Mesh file : mesh_vlm.dat");

    //Geometry changes for RANS & Structure
    // WING_LEFT.mirror_body(); 
    WING_RIGHT.change_all_distributions("partial", "cartesian");
    std::vector<std::vector<std::vector<std::vector<double>>>> surfaces_R = WING_RIGHT.get_paired_body_surfaces();    
    std::vector<std::vector<std::vector<std::vector<double>>>> surfaces_L = WING_LEFT.get_paired_body_surfaces();

    //RANS/Euler-----------
    
    std::vector<double> disc{100,150,200};

    std::cout<<"Surface size"<<std::endl;
    std::cout<<"start mesh rans"<<std::endl;

    for (int i=0; i < surfaces_R.size(); i++){
        file_name = {profile_name[i]+"_coarse.msh",profile_name[i]+"_normal.msh",profile_name[i]+"_fine.msh"};
        for (int j=0; j<3; j++){
            generer(surfaces_R[i][0], surfaces_R[i][1], surfaces_R[i][2], surfaces_R[i][4], surfaces_R[i][5], disc[j], file_name[j], viscous);
        }
    }  
    std::cout<<"end mesh rans"<<std::endl;
    gui.msg.push("[GEOM] Mesh for Euler/RANS solver generated");
    for(int i=0;i<file_name.size();i++){
        gui.msg.push("[GEOM] Mesh file :" + file_name[i]);
    }

    // Structure-------------
    std::vector<std::tuple<int,std::vector<double>,std::vector<double>,std::vector<double>>> element = maillage_structure(surfaces_R, settings.E, settings.G, settings.P_beam);
    gui.msg.push("[GEOM] Mesh for structural solver generated");
    gui.msg.push("[GEOM] Mesh file : Point_maillage_structure.txt");
}

void Geom::fill_database(database::table &table, std::vector<double> alphas){
        std::cout<<"Nom du profil database : ";
        std::cout<<profile_name[0]<<std::endl;


        for (int i=0; i<profile_name.size(); i++) {
            table.airfoils[profile_name[i]];                            // Créer le airfoil
            table.airfoils[profile_name[i]].alpha = {alphas};      //Remplir le champs alpha , 

            table.sectionAirfoils[i];                                   //WR
            table.sectionAirfoils[i] = {profile_name[i], profile_name[i]};   

            table.sectionSpanLocs[i];
            table.sectionSpanLocs[i] = {0.0,1.0};                       //doit aller de 0 à 1

            table.sectionAirfoils[i+1];                                 //WL
            table.sectionAirfoils[i+1] = {profile_name[i], profile_name[i]};   

            table.sectionSpanLocs[i+1];
            table.sectionSpanLocs[i+1] = {0.0,1.0}; 
        }
}

void Settings::import_config_file(tiny::config &io) {
	cr = io.get<double>("Geom-Wing", "Chord_root");
	ct = io.get<double>("Geom-Wing", "Chord_tip");
	envergure = io.get<double>("Geom-Wing", "Span");
	twist = io.get<double>("Geom-Wing", "Twist_angle");
	fleche = io.get<double>("Geom-Wing", "Sweep_angle");
	dihedre = io.get<double>("Geom-Wing", "Diherdral_angle");
    P_beam = io.get<double>("Geom-Wing", "Beam_position");
	//P_aile = io.get<double>("Geom-Wing", "Wing_position");
	Winglet = io.get<int>("Geom-Wing", "Winglet");
	S_type = io.get<int>("Geom-Wing", "Airfoil_type");

    if (solver_type() == "NACA") {
        c = io.get<double>("Geom-NACA", "c");
        p = io.get<double>("Geom-NACA", "p");
        t = io.get<double>("Geom-NACA", "t");
    } else if (solver_type() == "CST") {
        z_te = io.get<double>("Geom-CST", "Trailing_edge_dist");
        r_le = io.get<double>("Geom-CST", "Leading_edge_r");
        Beta = io.get<double>("Geom-CST", "Trailing_edge_angle");
    }

    E = io.get<double>("Geom-Meca", "E");
    G = io.get<double>("Geom-Meca", "G");
    mesh_VLM = io.get<bool>("Geom-Mesh","mesh_VLM");
    mesh_struct = io.get<bool>("Geom-Mesh","mesh_struct");
}

static const std::string bool_to_string(const bool b) {
    return b ? "true" : "false";
}

void Settings::export_config_file(tiny::config &io) {

	io.config["Geom-Wing"]["Chord_root"] = std::to_string(cr);
	io.config["Geom-Wing"]["Chord_tip"] = std::to_string(ct);
	io.config["Geom-Wing"]["Span"] = std::to_string(envergure);
	io.config["Geom-Wing"]["Twist_angle"] = std::to_string(twist);
	io.config["Geom-Wing"]["Sweep_angle"] = std::to_string(fleche);
	io.config["Geom-Wing"]["Diherdral_angle"] = std::to_string(dihedre);
	io.config["Geom-Wing"]["Beam_position"] = std::to_string(P_beam);
	//io.config["Geom-Wing"]["Wing_position"] = std::to_string(P_aile);
	io.config["Geom-Wing"]["Winglet"] = std::to_string(Winglet);
	io.config["Geom-Wing"]["Airfoil_type"] = std::to_string(S_type);
    
    if (solver_type() == "NACA") {
        io.config["Geom-NACA"]["m"] = std::to_string(c);
        io.config["Geom-NACA"]["p"] = std::to_string(p);
        io.config["Geom-NACA"]["t"] = std::to_string(t);
    } else if (solver_type() == "CST") {
        io.config["Geom-CST"]["Trailing_edge_dist"] = std::to_string(z_te);
        io.config["Geom-CST"]["Leading_edge_r"] = std::to_string(r_le);
        io.config["Geom-CST"]["Trailing_edge_angle"] = std::to_string(Beta);
    }

    io.config["Geom-Meca"]["E"] = std::to_string(E);
    io.config["Geom-Meca"]["G"] = std::to_string(G);
    io.config["Geom-Mesh"]["mesh_VLM"] = bool_to_string (mesh_VLM);
    io.config["Geom-Mesh"]["mesh_VLM"] = bool_to_string (mesh_struct);

}

