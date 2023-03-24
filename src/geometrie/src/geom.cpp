
#include "geometrie/geom.hpp"

using namespace geom;

Geom::Geom(GUIHandler &gui): gui(gui) {}

void Geom::Geom_gen() {
    std::string body_type = "General";
    Body WING_RIGHT(body_type);
    Body WING_LEFT(body_type);
    if(settings.S_type){
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
        {0, 2},                         //??, x)le                    
        {0.05, 0.1, 0.4},               //??, z_n            
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
        {0, 2},                         //??, x)le                    
        {0.05, 0.1, 0.4},               //??, z_n            
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
        {0, 2},                         //??, x)le 
        {0.05, 0.1, 0.4},               //??, z_n 
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
        {0, 2},                         //??, x)le 
        {0.05, 0.1, 0.4},               //??, z_n 
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
