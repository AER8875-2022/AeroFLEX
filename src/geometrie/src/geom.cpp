#include "geometrie/geom.hpp"
#include <cmath>
#include <iostream>
#include <string.h>

#ifndef M_PI
#define M_PI 3.141592653589793115997963468544185161590576171875
#endif

using namespace geom;

Geom::Geom(GUIHandler &gui): gui(gui) {}

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
    std::string body_type = "General";
    Body WING_RIGHT(body_type);
    Body WING_LEFT(body_type);
    int nb_profils = 1;
    //double perc_span = 1.0;
    for (int k=0; k<nb_profils; k++){
        if(settings.S_type == 1){
            //WR
            //#1
            WING_RIGHT.add_wing_surface_cst(
            100,                            
            5,
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
            5,
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
            profile_name = "CST";
        } else {
            //WR
            WING_RIGHT.add_wing_surface_naca(
            100,                            
            5,                              
            settings.envergure,             
            0, 
            {settings.cr,settings.ct}, 
            {settings.m},                   
            {settings.p},                        
            {settings.t},                     
            {0.0,xle},                                           
            {0.0,zn},                          
            {0.0,d_alpha}                      
            );

            //WL
            WING_LEFT.add_wing_surface_naca(
            100,                            
            5,                              
            settings.envergure,             
            0, 
            {settings.cr,settings.ct}, 
            {settings.m},                   
            {settings.p},                        
            {settings.t},                     
            {0.0,xle},                                           
            {0.0,zn},                          
            {0.0,d_alpha}                     
            );
            profile_name = "NACA";
        }
    }

    WING_LEFT.mirror_body(); 
    WING_RIGHT.change_all_distributions("partial", "cartesian");
    WR_surfaces = WING_RIGHT.get_paired_body_surfaces();    
    std::vector<std::vector<std::vector<std::vector<double>>>> WL_surfaces = WING_LEFT.get_paired_body_surfaces();
    gui.msg.push("[GEOM] Geometry generated");


    // Structure
    std::vector<std::tuple<int,std::vector<double>,std::vector<double>,std::vector<double>>> element = maillage_structure(WING_RIGHT, settings.E, settings.G);
    gui.msg.push("[GEOM] Mesh for structure solver generated");
}

//RANS
void Geom::Geom_mesh(bool viscous) {
    std::vector<double> disc{100,150,200};
    std::vector<std::vector<std::vector<std::vector<double>>>> surfaces = WR_surfaces;
    file_name = {profile_name+"_coarse.msh",profile_name+"_normal.msh",profile_name+"_fine.msh"};

    //std::cout<<"Surface size"<<std::endl;
    std::cout<<"start mesh rans"<<std::endl;

    for (int i=0; i < surfaces.size(); i++){
        gui.msg.push("[GEOM] test");
        for (int j=0; j<3; j++){
            generer(surfaces[i][0], surfaces[i][1], surfaces[i][2], surfaces[i][4], surfaces[i][5], disc[j], file_name[j], viscous);
        }
    }  
    std::cout<<"end mesh rans"<<std::endl;
    gui.msg.push("[GEOM] Mesh for Euler/RANS solver generated");
}

void Geom::fill_database(database::table &table){
        std::cout<<"Nom du profil database :";
        std::cout<<profile_name<<std::endl;
        table.airfoils[profile_name]; // Créer le airfoil
        table.airfoils[profile_name].alpha = {0.0,5.0,10.0};    //Remplir le champs alpha 

        table.sectionAirfoils[0];       //0 aile droit, 1 aile gauche
        table.sectionAirfoils[0] = {profile_name, profile_name};   

        table.sectionSpanLocs[0];
        table.sectionSpanLocs[0] = {0.0,1.0}; //doit aller de 0 à 1
}

void Settings::import_config_file(tiny::config &io) {
	cr = io.get<double>("Geom-Wing", "Chord_root");
	ct = io.get<double>("Geom-Wing", "Chord_tip");
	envergure = io.get<double>("Geom-Wing", "Span");
	twist = io.get<double>("Geom-Wing", "Twist_angle");
	fleche = io.get<double>("Geom-Wing", "Sweep_angle");
	dihedre = io.get<double>("Geom-Wing", "Diherdral_angle");
    P_beam = io.get<double>("Geom-Wing", "Beam_position");
	P_aile = io.get<double>("Geom-Wing", "Wing_position");
	Winglet = io.get<int>("Geom-Wing", "Winglet");
	S_type = io.get<int>("Geom-Wing", "Airfoil_type");

    if (solver_type() == "NACA") {
        m = io.get<double>("Geom-NACA", "m");
        p = io.get<double>("Geom-NACA", "p");
        t = io.get<double>("Geom-NACA", "t");
    } else if (solver_type() == "CST") {
        z_te = io.get<double>("Geom-CST", "Trailing_edge_dist");
        r_le = io.get<double>("Geom-CST", "Leading_edge_r");
        Beta = io.get<double>("Geom-CST", "Trailing_edge_angle");
    }

    E = io.get<double>("Geom-Meca", "E");
    G = io.get<double>("Geom-Meca", "G");
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
	io.config["Geom-Wing"]["Wing_position"] = std::to_string(P_aile);
	io.config["Geom-Wing"]["Winglet"] = std::to_string(Winglet);
	io.config["Geom-Wing"]["Airfoil_type"] = std::to_string(S_type);
    
    if (solver_type() == "NACA") {
        io.config["Geom-NACA"]["m"] = std::to_string(m);
        io.config["Geom-NACA"]["p"] = std::to_string(p);
        io.config["Geom-NACA"]["t"] = std::to_string(t);
    } else if (solver_type() == "CST") {
        io.config["Geom-CST"]["Trailing_edge_dist"] = std::to_string(z_te);
        io.config["Geom-CST"]["Leading_edge_r"] = std::to_string(r_le);
        io.config["Geom-CST"]["Trailing_edge_angle"] = std::to_string(Beta);
    }

    io.config["Geom-Meca"]["E"] = std::to_string(E);
    io.config["Geom-Meca"]["G"] = std::to_string(G);
}

