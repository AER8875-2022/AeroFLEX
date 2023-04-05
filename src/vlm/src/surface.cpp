
#include "vlm/surface.hpp"
#include "iostream"
#include <cmath>

using namespace vlm;
using namespace surface;
using namespace Eigen;

wingStation::wingStation(const int globalIndex,
                         const std::vector<int> &vortexIDs)
    : globalIndex(globalIndex), vortexIDs(vortexIDs), area(0.0), cl(0.0),
      cm(Vector3d::Zero()) {}

void wingStation::initialize(const std::vector<Vector3d> &nodes,
                             std::vector<element::vortexRing> &vortices,
                             const input::simParam &sim) {
  for (auto &vortexID : vortexIDs) {
    vortices[vortexID].initialize(nodes, sim);
  }
  computeArea(vortices);
  computeChordLength(nodes, vortices);
  // Computing the undeformed spanwise location for interpolation
  spanLoc = forceActingPoint(nodes, vortices)(1);
  // Setting the local aoa to the geometric aoa
  local_aoa = sim.aoa;
}

void wingStation::generateWake(const double wakeLength,
                               const std::vector<Vector3d> &nodes,
                               const std::vector<element::vortexRing> &vortices,
                               const input::simParam &sim,
                               std::vector<Vector3d> &wakeNodes,
                               std::vector<element::vortexRing> &wakePanels) {
  auto &trailPanel = vortices[vortexIDs.back()];
  // Adding trailing edge nodes
  auto p1 = nodes[trailPanel.get_nodeIDs()[3]];
  auto p2 = nodes[trailPanel.get_nodeIDs()[2]];
  wakeNodes.push_back(p1);
  wakeNodes.push_back(p2);
  // Generating new wake nodes
  auto &wakeOrientation = sim.freeStream().normalized();
  wakeNodes.push_back(p2 + wakeLength * wakeOrientation);
  wakeNodes.push_back(p1 + wakeLength * wakeOrientation);
  // Generating new wake panel
  int globalIndex = trailPanel.get_globalIndex();
  int id1 = wakeNodes.size() - 4;
  int id2 = wakeNodes.size() - 3;
  int id3 = wakeNodes.size() - 2;
  int id4 = wakeNodes.size() - 1;
  auto wakePanel = element::vortexRing(globalIndex, {id1, id2, id3, id4}, 1.0);
  wakePanel.initialize(wakeNodes, sim);
  wakePanels.push_back(wakePanel);
}

void wingStation::updateGeometry(const std::vector<Vector3d> &nodes,
                                 std::vector<element::vortexRing> &vortices) {
  for (auto &vortexID : vortexIDs) {
    vortices[vortexID].updateGeometry(nodes);
  }
  computeArea(vortices);
  computeChordLength(nodes, vortices);
}

Vector3d wingStation::forceActingPoint(
    const std::vector<Vector3d> &nodes,
    const std::vector<element::vortexRing> &vortices) const {
  return (vortices[vortexIDs.front()].forceActingPoint(nodes));
}

void wingStation::updateLocalAoa(const double dalpha,
                                 const std::vector<Vector3d> &nodes,
                                 std::vector<element::vortexRing> &vortices) {
  local_aoa += dalpha;
  for (auto &vortexID : vortexIDs) {
    vortices[vortexID].local_aoa += dalpha;
    vortices[vortexID].computeCollocationPoint(nodes);
  }
}

double wingStation::get_globalIndex() const { return globalIndex; }

double wingStation::get_area() const { return area; }

double wingStation::get_chord() const { return chord; }

double wingStation::get_spanLoc() const { return spanLoc; }

std::vector<int> wingStation::get_vortexIDs() const { return vortexIDs; }

double wingStation::get_cl() const { return cl; }

double wingStation::get_cd() const { return cd; }

Vector3d wingStation::get_cm() const { return cm; }

void wingStation::computeArea(
    const std::vector<element::vortexRing> &vortices) {
  area = 0.0;
  for (auto &vortexID : vortexIDs) {
    area += vortices[vortexID].get_area();
  }
}

void wingStation::computeChordLength(
    const std::vector<Vector3d> &nodes,
    const std::vector<element::vortexRing> &vortices) {
  auto &leadingEdge = vortices[vortexIDs.front()];
  auto &trailingEdge = vortices[vortexIDs.back()];

  double d1 = nodes[trailingEdge.get_nodeIDs()[3]](0) -
              nodes[leadingEdge.get_nodeIDs()[0]](0);
  double d2 = nodes[trailingEdge.get_nodeIDs()[2]](0) -
              nodes[leadingEdge.get_nodeIDs()[1]](0);
  chord = 0.5 * (d1 + d2);
}

void wingStation::computeForces(const input::simParam &sim,
                                const std::vector<Vector3d> &nodes,
                                std::vector<element::vortexRing> &vortices) {
  Vector3d Qinf = sim.freeStream(local_aoa);
  Vector3d liftAxis = Qinf.cross(Vector3d::UnitY()).normalized();
  double lift = 0.0;
  Vector3d moment = Vector3d::Zero();
  double previousGamma = 0.0;
  for (auto &vortexID : vortexIDs) {
    auto &vortex = vortices[vortexID];

    // Computing distances
    Vector3d dl = vortex.leadingEdgeDl(nodes);
    Vector3d lever = sim.origin - vortex.forceActingPoint(nodes);

    // Vectorial force
    Vector3d force = Qinf.cross(dl);
    force *= sim.rho * (vortex.get_gamma() - previousGamma);

    // Computing local force and moment
    double localForce = force.dot(liftAxis);
    Vector3d localMoment = force.cross(lever);

    // Saving forces of the panel
    vortex.cl =
        localForce / (0.5 * sim.rho * Qinf.norm() * Qinf.norm() * sim.sref);
    vortex.cm = localMoment / (0.5 * sim.rho * Qinf.norm() * Qinf.norm() *
                               sim.cref * sim.sref);

    // Incrementing wing station forces
    lift += localForce;
    moment += localMoment;

    // Setting gamma reference for next panel
    previousGamma = vortex.get_gamma();
  }
  // Computing coefficients
  cl = lift / (0.5 * sim.rho * Qinf.norm() * Qinf.norm() * area);
  cm = moment / (0.5 * sim.rho * Qinf.norm() * Qinf.norm() * area * chord);
}

// ----------------------------

wing::wing(const int globalIndex, const std::vector<int> &stationIDs)
    : globalIndex(globalIndex), stationIDs(stationIDs), area(0.0), cl(0.0),
      cm(Vector3d::Zero()) {}

void wing::initialize(const std::vector<Vector3d> &nodes,
                      std::vector<wingStation> &stations,
                      std::vector<element::vortexRing> &vortices,
                      const input::simParam &sim) {
  for (auto &stationID : stationIDs) {
    stations[stationID].initialize(nodes, vortices, sim);
  }
  computeArea(stations);
  // Computing span of the wing
  auto &rootStation = stations[stationIDs.front()];
  auto &tipStation = stations[stationIDs.back()];
  span = tipStation.forceActingPoint(nodes, vortices)(1) -
         rootStation.forceActingPoint(nodes, vortices)(1);
}

void wing::updateGeometry(const std::vector<Vector3d> &nodes,
                          std::vector<wingStation> &stations,
                          std::vector<element::vortexRing> &vortices) {
  for (auto &stationID : stationIDs) {
    stations[stationID].updateGeometry(nodes, vortices);
  }
  computeArea(stations);
}

double wing::get_globalIndex() const { return globalIndex; }

double wing::get_area() const { return area; }

double wing::get_span() const { return span; }

std::vector<int> wing::get_stationIDs() const { return stationIDs; }

double wing::get_cl() const { return cl; }

double wing::get_cd() const { return cd; }

Vector3d wing::get_cm() const { return cm; }

void wing::computeArea(const std::vector<wingStation> &stations) {
  area = 0.0;
  for (auto &stationID : stationIDs) {
    area += stations[stationID].get_area();
  }
}

void wing::computeForces(const input::simParam &sim,
                         const std::vector<surface::wingStation> &stations) {
  cl = 0.0;
  cd = 0.0;
  cm = Vector3d::Zero();
  for (auto &stationID : stationIDs) {
    auto &station = stations[stationID];
    cl += station.get_cl() * station.get_area() / area;
    cd += station.get_cd() * station.get_area() / area;
    cm += station.get_cm() * station.get_area() / area * station.get_chord() /
          sim.cref;
  }
}

// ----------------------------

patch::patch(const int globalIndex, const std::vector<int> &doubletIDs)
    : globalIndex(globalIndex), doubletIDs(doubletIDs), area(0.0) {}

void patch::initialize(const std::vector<Vector3d> &nodes,
                       std::vector<element::doubletPanel> &doublets,
                       const input::simParam &sim) {
  for (auto &doubletID : doubletIDs) {
    doublets[doubletID].initialize(nodes, sim);
  }
  computeArea(doublets);
}

void patch::updateGeometry(const std::vector<Vector3d> &nodes,
                           std::vector<element::doubletPanel> &doublets) {
  for (auto &doubletID : doubletIDs) {
    doublets[doubletID].updateGeometry(nodes);
  }
  computeArea(doublets);
}

void patch::ScanNeighbor(std::vector<element::doubletPanel> &doublets) {
  //for the target panel
  for(auto &doub_target : doublets) {
    auto nodes_target = doub_target.get_nodeIDs();
    //std::cout<< "new panel "<< doub_target.get_globalIndex() << " : " << std::endl;
    for (auto &doub_compared : doublets){
      if (doub_target.get_globalIndex() != doub_compared.get_globalIndex()){
        auto nodes_compared = doub_compared.get_nodeIDs();
        for (size_t i=0; i<nodes_target.size()-1; i++){
          //search for the first point of each edge
          auto it1 = std::find(nodes_compared.begin(), nodes_compared.end(), nodes_target[i]);
          //search for the second point of each value
          auto it2 = std::find(nodes_compared.begin(), nodes_compared.end(), nodes_target[i+1]);
          //if value found search for the second
          if( it1 != nodes_compared.end() && it2 != nodes_compared.end()) {
            doub_target.NeighborPanel_IDs[i] = doub_compared.get_globalIndex();
            //std::cout<< doub_target.NeighborPanel_IDs[i]<< std::endl;
          } else {

          } //******NEED TO CORRECT FOR WHEN POINTS ARE COINCIDED******
        }
        //for the last point  
        auto it1 = std::find(nodes_compared.begin(), nodes_compared.end(), nodes_target.back());
          //search for the second point of each value
        auto it2 = std::find(nodes_compared.begin(), nodes_compared.end(), nodes_target.front());
        if( it1 != nodes_compared.end() && it2 != nodes_compared.end()) {
           doub_target.NeighborPanel_IDs[nodes_target.size()-1] = doub_compared.get_globalIndex();
          //std::cout<< doub_target.NeighborPanel_IDs[nodes_target.size()-1]<< std::endl;
        } else {
            
        }
      }
    }
  }
}

//*******not used removed during clean up *********
void patch::Storing_nondirectPanel(std::vector<element::doubletPanel> &doublets) {
  //int count=0;
  for (auto &doub : doublets){
    auto neighbor = doub.get_neighbor();
    for (size_t i=0; i<neighbor.size(); i++){
    //auto it = std::find(neighbor.begin(), neighbor.end(), -1);

      //methode pas efficace
      if(neighbor[i] == -1){
        //doub.nondirectNeighbor_IDs[i]= doub[neighbor[i-2]].NeighborPanel_IDs[i];
      }
      else {
        doub.nondirectNeighbor_IDs[i]= -1 ;
      }
    }
    //count++;
  }
}

void patch::computePressure(const input::simParam &sim,
               std::vector<element::doubletPanel> &doublets) {
  cp = 0.0; //remove in clean up
  const auto &vinf = sim.vinf;
  const Vector3d freestream = sim.freeStream();
  for (auto &doublet : doublets) {
    Vector3d sum = Vector3d::Zero();
    auto neighbor = doublet.get_neighbor();
    auto edge_center = doublet.get_edge_center();
    auto edge = doublet.get_edges(); // getting the length of the edge from the mesh (not the co-planair panel) [possible mistake]
    auto area = doublet.get_area();
    //Computing Velocity using Green-Gauss
    //std::cout << "Panel number  " << doublet.get_globalIndex() << std::endl;
    for (size_t i=0; i<neighbor.size(); i++){
      if (neighbor[i]!=-1){
        double f = (doublets[neighbor[i]].get_center() - edge_center[i]).norm() / (doublets[neighbor[i]].get_center() - doublet.get_center()).norm(); 
        double phi_f = (doublet.mu * f) + (doublets[neighbor[i]].mu * (1 - f));
        //double phi_f = (12 * f) + (7 * (1 - f)); //for debugging
        //doublet.get_normal()
        
        Vector3d local_panel_normal = {0, 0, 1};
        Vector3d edge_vector_local = doublet.ProjectingToLocal(edge[i].get_dl());
        //edge_vector_local[2] = 0;
        Vector3d segment_normal = (local_panel_normal).cross(edge_vector_local);
        
        /*
        std::cout << "segment_normal :" << std::endl;
            for (int w=0 ; w<3; w++){
            std::cout << segment_normal[w] << std::endl;
          } 
          std::cout << "phi_f :" << std::endl;
          std::cout << phi_f << std::endl;
        */

        if (i==0){
          doublet.segment_normal = segment_normal.normalized();
        }
        /*
        std::cout << "edge numero " << i << " :  " << std::endl;
        std::cout << "center of the edge " << ":  " << std::endl;
        std::cout << edge_center[i] << std::endl; 
        std::cout << "weighing geometric factor (f) : ";
        std::cout << f <<std::endl;
        std::cout << "Value of MU taking account the geometry (Phi_f) : ";
        std::cout << phi_f <<std::endl;
        */
        sum += phi_f * segment_normal.normalized() * edge[i].get_dl().norm();
      }
    }
      
    //La mauvaise aire est utilisé ici. Il faut plutôt prendre l'aire du panneau co-planaire
    Vector3d gradMU = (1 / area) * sum;
    Vector3d localstream = doublet.ProjectingToLocal(freestream); //
    localstream[2]=0;
    doublet.localstream = localstream; //for troubleshooting
    gradMU+=localstream;

    /*
    std::cout << "velocity :" << std::endl;
      for (int w=0 ; w<3; w++){
        std::cout << gradMU[w] << std::endl;
      } 
    */
    double V_k = sqrt(gradMU[0]*gradMU[0] + gradMU[1]*gradMU[1]);
    if (doublet.get_globalIndex()>=110 && doublet.get_globalIndex()<=119){

      //std::cout<< doublet.get_globalIndex() << std::endl;
      std::cout<< V_k/vinf << std::endl;
      //std::cout << doublet.get_center()[0] << std::endl;
             
    }

    doublet.storing_velocity(gradMU, V_k/vinf);
    //std::cout << "vitesse K :  " << V_k << std::endl; 
    doublet.cp = 1 - (V_k*V_k)/(vinf*vinf);
    //std::cout << "coefficient de pression (Cp)  : ";
    //std::cout << doublet.cp << std::endl;
  }      
}

double patch::get_globalIndex() const { return globalIndex; }

double patch::get_area() const { return area; }

double patch::get_cp() const { return cp; }

std::vector<int> patch::get_doubletIDs() const { return doubletIDs; }

void patch::computeArea(const std::vector<element::doubletPanel> &doublets) {
  area = 0.0;
  for (auto &doubletID : doubletIDs) {
    area += doublets[doubletID].get_area();
  }
}
