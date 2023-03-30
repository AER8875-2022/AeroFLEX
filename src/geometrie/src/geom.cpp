
#include "geometrie/geom.hpp"
#include "geometrie/euler.hpp"
#include "geometrie/geometry.hpp"
#include "geometrie/structure.hpp"

using namespace geom;

Geom::Geom(GUIHandler &gui): gui(gui) {}

void Geom::Geom_gen() {
    std::string body_type = "General";
    Body WING_RIGHT(body_type);
    Body WING_LEFT(body_type);
    if(settings.S_type == 1){
        //WR
        WING_RIGHT.add_wing_surface_cst(
        100,                            
        5,
        settings.envergure,             
        0,                              
        {settings.cr,settings.ct},      
        {settings.z_te},                         
        {settings.r_le},                        
        {settings.Beta},                         
        {0.0},                         //??, x)le                    
        {0.0},               //??, z_n            
        {settings.twist}                         
        );
        //WL
        WING_LEFT.add_wing_surface_cst(
        100,                            
        5,
        settings.envergure,             
        0,                              
        {settings.cr,settings.ct},      
        {settings.z_te},                         
        {settings.r_le},                        
        {settings.Beta},                         
        {0.0},                         //??, x)le                    
        {0.0},               //??, z_n            
        {settings.twist}                         
        );

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
        {0.0},                         //??, x)le 
        {0.0},               //??, z_n 
        {settings.twist}                      
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
        {0.0},                         //??, x)le 
        {0.0},               //??, z_n 
        {settings.twist}                      
        );
    }

    WING_LEFT.mirror_body(); 
    WING_RIGHT.change_all_distributions("partial", "cartesian");
    std::vector<std::vector<std::vector<std::vector<double>>>> WR_surfaces = WING_RIGHT.get_paired_body_surfaces();    
    std::vector<std::vector<std::vector<std::vector<double>>>> WL_surfaces = WING_LEFT.get_paired_body_surfaces();
    
    //Eventuellement changer de fonction

    //RANS
    std::vector<double> disc{100,150,200};
    //cout<<disc[0]<<endl;
    std::vector<std::string> file_name{"Airfoil_coarse.msh","Airfoil_normal.msh","Airfoil_fine.msh"};
    //std::cout<<file_name[0]<<endl;
    std::vector<std::vector<std::vector<std::vector<double>>>> surfaces = WR_surfaces;
    bool RANS = false; // True = solver RANS, False = solver Euler
    for (int i=0; i < surfaces.size(); i++){
        for (int j=0; j<3; j++){
            generer(surfaces[i][0], surfaces[i][1], surfaces[i][2], surfaces[i][4], surfaces[i][5], disc[j], file_name[j], RANS);
        }
    }  

    // Structure
    std::vector<std::tuple<int,std::vector<double>,std::vector<double>,std::vector<double>>> element = maillage_structure(WING_RIGHT);
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

