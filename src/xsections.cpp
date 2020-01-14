 /******************************************************************************
 *    This program is free software: you can redistribute it and/or modify     *
 *   it under the terms of the GNU General Public License as published by      *
 *   the Free Software Foundation, either version 3 of the License, or         *
 *   (at your option) any later version.                                       *
 *                                                                             *
 *   This program is distributed in the hope that it will be useful,           *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *   GNU General Public License for more details.                              *
 *                                                                             *
 *   You should have received a copy of the GNU General Public License         *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.     *
 *                                                                             *
 *   Authors:                                                                  *
 *      Carlos Arguelles (University of Wisconsin Madison)                     *
 *         carguelles@icecube.wisc.edu                                         *
 *      Jordi Salvado (University of Wisconsin Madison)                        *
 *         jsalvado@icecube.wisc.edu                                           *
 *      Christopher Weaver (University of Wisconsin Madison)                   *
 *         chris.weaver@icecube.wisc.edu                                       *
 ******************************************************************************/

#include <nuSQuIDS/xsections.h>
#include <iostream>
#include <fstream>

namespace{
struct H5File{
  hid_t id;
  H5File(hid_t id):id(id){}
  H5File(const H5File&)=delete;
  H5File(H5File&& h):id(h.id){ h.id=0; }
  H5File& operator=(const H5File&)=delete;
  H5File& operator=(H5File&& h){ std::swap(id,h.id); return *this;}
  ~H5File(){ H5Fclose(id); }
  operator hid_t() const{ return(id); }
};
}

namespace nusquids{

double NeutrinoDISCrossSectionsFromTables::LinInter(double x,double xM, double xP,double yM,double yP) const{
  return yM + (yP-yM)*(x-xM)/(xP-xM);
}

double NeutrinoDISCrossSectionsFromTables::TotalCrossSection(double Enu, NeutrinoFlavor flavor,
                           NeutrinoType neutype, Current current, Target target) const{
  // we assume that sterile neutrinos are truly sterile
  if (not (flavor == electron or flavor == muon or flavor == tau))
    return 0.0;

  if (Enu > Emax)
    throw std::runtime_error("NeutrinoCrossSections::TotalCrossSection: Only DIS cross sections are included in a limited range. Interpolation re\
quested below "+std::to_string(Emin/GeV)+" GeV or above "+std::to_string(Emax/GeV)+" GeV. E_nu = " + std::to_string(Enu/GeV) + " [GeV].");

  if (Enu < 10*GeV and Emin > 10*GeV) {
    std::clog << "NeutrinoCrossSections::TotalCrossSection: Neglecting the neutrino cross section below 10 GeV." << std::endl;
    return std::numeric_limits<double>::min();
  } else if (Enu < Emin) {
    // use approximate linear scaling below Emin GeV which is a good approximation for DIS up to 10 GeV
    return TotalCrossSection(Emin, flavor, neutype, current, target)*(Enu/Emin);
  }

  // convert to GeV
  Enu /= GeV;

  double logE = log(Enu);
  double dlogE = logE_data_range[1]-logE_data_range[0];
  size_t idx = static_cast<size_t>((logE-logE_data_range[0])/dlogE);
  if(idx==div)
    idx--;
  const marray<double,4>& sigma=(current==CC ? s_CC_data : s_NC_data);
  return(LinInter(logE,logE_data_range[idx],logE_data_range[idx+1],sigma[target][neutype][flavor][idx],sigma[target][neutype][flavor][idx+1]));
}

double NeutrinoDISCrossSectionsFromTables::SingleDifferentialCrossSection(double E1, double E2, NeutrinoFlavor flavor, NeutrinoType neutype, Current current, Target target) const{
  // we assume that sterile neutrinos are trully sterile
  if (not (flavor == electron or flavor == muon or flavor == tau))
    return 0.0;

  if (E1 > Emax)
    throw std::runtime_error("NeutrinoCrossSections::SingleDifferentialCrossSection: Only DIS cross sections are included. Interpolation re\
quested below "+std::to_string(Emin/GeV)+" GeV or above "+std::to_string(Emax/GeV)+" GeV. E_nu = " + std::to_string(E1/GeV) + " [GeV].");
  if (E1 < Emin)
    return std::numeric_limits<double>::min();

  // convert to GeV
  E1 /= GeV;
  E2 /= GeV;

  double logE1 = log(E1);
  double logE2 = log(E2);
  double dlogE = logE_data_range[1]-logE_data_range[0];

  size_t loge_M1 = static_cast<size_t>((logE1-logE_data_range[0])/dlogE);
  size_t loge_M2 = static_cast<size_t>((logE2-logE_data_range[0])/dlogE);

  if ((loge_M2 > div-1) or (loge_M1 > div-1) or (E2 >= E1))
    return 0.0;

  //std::cout << E1 << " " << E2 << " " << loge_M1 << " " << loge_M2 << " " << div << std::endl;
  double phiMM,phiMP,phiPM,phiPP;
  if (current == CC){
    phiMM = dsde_CC_data[target][neutype][flavor][loge_M1][loge_M2];
    phiMP = dsde_CC_data[target][neutype][flavor][loge_M1][loge_M2+1];
    if ( loge_M1 == div-1 ){
      // we are at the boundary, cannot bilinearly interpolate
      return LinInter(logE2,logE_data_range[loge_M2],logE_data_range[loge_M2+1],
          phiMM,phiMP);
    }
    phiPM = dsde_CC_data[target][neutype][flavor][loge_M1+1][loge_M2];
    phiPP = dsde_CC_data[target][neutype][flavor][loge_M1+1][loge_M2+1];
  } else if (current == NC){
    phiMM = dsde_NC_data[target][neutype][flavor][loge_M1][loge_M2];
    phiMP = dsde_NC_data[target][neutype][flavor][loge_M1][loge_M2+1];
    if ( loge_M1 == div-1 ){
      // we are at the boundary, cannot bilinearly interpolate
      return LinInter(logE2,logE_data_range[loge_M2],logE_data_range[loge_M2+1],
          phiMM,phiMP);
    }
    phiPM = dsde_NC_data[target][neutype][flavor][loge_M1+1][loge_M2];
    phiPP = dsde_NC_data[target][neutype][flavor][loge_M1+1][loge_M2+1];
  } else
    throw std::runtime_error("nuSQUIDS::XSECTIONS::ERROR::Current type unkwown.");

  return LinInter(logE1,logE_data_range[loge_M1],logE_data_range[loge_M1+1],
           LinInter(logE2,logE_data_range[loge_M2],logE_data_range[loge_M2+1],phiMM,phiMP),
           LinInter(logE2,logE_data_range[loge_M2],logE_data_range[loge_M2+1],phiPM,phiPP));
}

void NeutrinoDISCrossSectionsFromTables::ReadText(std::string root){

       // Define proton and neutron cross section filenames
       std::string filename_p_dsde_CC = root+"p_dsde_CC.dat";
       std::string filename_p_dsde_NC = root+"p_dsde_NC.dat";
       std::string filename_n_dsde_CC = root+"n_dsde_CC.dat";
       std::string filename_n_dsde_NC = root+"n_dsde_NC.dat";
       std::string filename_p_sigma_CC = root+"p_sigma_CC.dat";
       std::string filename_p_sigma_NC = root+"p_sigma_NC.dat";
       std::string filename_n_sigma_CC = root+"n_sigma_CC.dat";
       std::string filename_n_sigma_NC = root+"n_sigma_NC.dat";

       // Define isoscalar cross section file names
       std::string filename_dsde_CC = root+"dsde_CC.dat";
       std::string filename_dsde_NC = root+"dsde_NC.dat";
       std::string filename_sigma_CC = root+"sigma_CC.dat";
       std::string filename_sigma_NC = root+"sigma_NC.dat";

       // Check if proton and neutron cross section files exist. If yes, use p/n by default
       if(
          fexists(filename_p_dsde_CC) and
          fexists(filename_p_dsde_NC) and
          fexists(filename_n_dsde_CC) and
          fexists(filename_n_dsde_NC) and
          fexists(filename_p_sigma_CC) and
          fexists(filename_p_sigma_NC) and
          fexists(filename_n_sigma_CC) and
          fexists(filename_n_sigma_NC)
         )
       {
          marray<double,2> p_dsde_CC_raw_data = quickread(filename_p_dsde_CC);
          marray<double,2> p_dsde_NC_raw_data = quickread(filename_p_dsde_NC);
          marray<double,2> n_dsde_CC_raw_data = quickread(filename_n_dsde_CC);
          marray<double,2> n_dsde_NC_raw_data = quickread(filename_n_dsde_NC);
          marray<double,2> p_sigma_CC_raw_data = quickread(filename_p_sigma_CC);
          marray<double,2> p_sigma_NC_raw_data = quickread(filename_p_sigma_NC);
          marray<double,2> n_sigma_CC_raw_data = quickread(filename_n_sigma_CC);
          marray<double,2> n_sigma_NC_raw_data = quickread(filename_n_sigma_NC);
          use_isoscalar = false;
       }

       // check if isoscalar cross section files exist
       else if (
                 fexists(filename_niso_dsde_CC) and
                 fexists(filename_niso_dsde_NC) and
                 fexists(filename_niso_sigma_CC) and
                 fexists(filename_niso_sigma_NC)
                )
       {
          marray<double,2> niso_dsde_CC_raw_data = quickread(filename_niso_dsde_CC);
          marray<double,2> niso_dsde_NC_raw_data = quickread(filename_niso_dsde_NC);
          marray<double,2> niso_sigma_CC_raw_data = quickread(filename_niso_sigma_CC);
          marray<double,2> niso_sigma_NC_raw_data = quickread(filename_niso_sigma_NC);
          use_isoscalar = true;
       }
       else {
         throw std::runtime_error("nuSQUIDS::XSECTIONS::ERROR::Cross section files not found.");
       }
       if(use_isoscalar){
          registered_targets_ = {Isoscalar};
/*          // read data tables
          marray<double,2> niso_dsde_CC_raw_data = quickread(filename_dsde_CC);
          marray<double,2> dsde_NC_raw_data = quickread(filename_dsde_NC);
          marray<double,2> sigma_CC_raw_data = quickread(filename_sigma_CC);
          marray<double,2> sigma_NC_raw_data = quickread(filename_sigma_NC);
*/
          // check table shapes and get the number of energy nodes
          unsigned int data_e_size = 0;
          if( niso_sigma_CC_raw_data.extent(0) == niso_sigma_NC_raw_data.extent(0) and
              niso_sigma_NC_raw_data.extent(1) == niso_sigma_NC_raw_data.extent(1) )
            data_e_size = niso_sigma_CC_raw_data.extent(0);
          else
            throw std::runtime_error("nuSQUIDS::xsections::init: Data tables not the same size.");

          // getting the raw data energy node values
          logE_data_range.resize(data_e_size);
          for( int ie = 0; ie < data_e_size; ie ++)
            logE_data_range[ie] = log(niso_sigma_CC_raw_data[ie][0]);

          Emin = niso_sigma_CC_raw_data[0][0]*GeV;
          Emax = niso_sigma_CC_raw_data[data_e_size-1][0]*GeV;
          div = data_e_size;

          // convert raw data tables into formatted marrays
          s_CC_data.resize(std::vector<size_t>{1,2,3,data_e_size});
          s_NC_data.resize(std::vector<size_t>{1,2,3,data_e_size});
          dsde_CC_data.resize(std::vector<size_t>{1,2,3,data_e_size,data_e_size});
          dsde_NC_data.resize(std::vector<size_t>{1,2,3,data_e_size,data_e_size});

          for(NeutrinoType neutype : {neutrino,antineutrino}){
            for(NeutrinoFlavor flavor : {electron,muon,tau}){
              for(unsigned int e1 = 0; e1 < data_e_size; e1++){
                  // FIXME
                s_CC_data[Isoscalar][neutype][flavor][e1] = niso_sigma_CC_raw_data[e1][1+2*((int)flavor)+(int)neutype];
                  s_NC_data[Isoscalar][neutype][flavor][e1] = niso_sigma_NC_raw_data[e1][1+2*((int)flavor)+(int)neutype];
                  for (unsigned int e2 = 0; e2 < data_e_size; e2++){
                    dsde_CC_data[Isoscalar][neutype][flavor][e1][e2] = niso_dsde_CC_raw_data[e1*data_e_size+e2][2+2*(static_cast<int>(flavor))+static_cast<int>(neutype)];
                    dsde_NC_data[Isoscalar][neutype][flavor][e1][e2] = niso_dsde_NC_raw_data[e1*data_e_size+e2][2+2*(static_cast<int>(flavor))+static_cast<int>(neutype)];
                  }
                }
              }
            }
      } else {
          registered_targets_ = {Proton, Neutron};
          unsigned int data_e_size = 0;
          if( p_sigma_CC_raw_data.extent(0) == p_sigma_NC_raw_data.extent(0) and
              p_sigma_CC_raw_data.extent(0) == p_sigma_NC_raw_data.extent(0) and 
              n_sigma_NC_raw_data.extent(1) == n_sigma_NC_raw_data.extent(1) and
              n_sigma_NC_raw_data.extent(1) == n_sigma_NC_raw_data.extent(1)
            ) {
                data_e_size = p_sigma_CC_raw_data.extent(0);
              }
          else
            throw std::runtime_error("nuSQUIDS::xsections::init: Data tables not the same size.");

          // getting the raw data energy node values
          logE_data_range.resize(data_e_size);
          for( int ie = 0; ie < data_e_size; ie ++)
            logE_data_range[ie] = log(p_sigma_CC_raw_data[ie][0]);

          // TODO find out if it's okay to use only proton info
          Emin = p_sigma_CC_raw_data[0][0]*GeV;
          Emax = p_sigma_CC_raw_data[data_e_size-1][0]*GeV;
          div = data_e_size;

          // convert raw data tables into formatted marrays
          s_CC_data.resize(std::vector<size_t>{2,2,3,data_e_size});
          s_NC_data.resize(std::vector<size_t>{2,2,3,data_e_size});

          dsde_CC_data.resize(std::vector<size_t>{2,2,3,data_e_size,data_e_size});
          dsde_NC_data.resize(std::vector<size_t>{2,2,3,data_e_size,data_e_size});
          for(Target target : {Proton, Neutron}){
            for(NeutrinoType neutype : {neutrino,antineutrino}){
              for(NeutrinoFlavor flavor : {electron,muon,tau}){
                for(unsigned int e1 = 0; e1 < data_e_size; e1++){
                  if(target == Proton){
                    s_CC_data[target][neutype][flavor][e1] = p_sigma_CC_raw_data[e1][1+2*((int)flavor)+(int)neutype];
                    s_NC_data[target][neutype][flavor][e1] = p_sigma_NC_raw_data[e1][1+2*((int)flavor)+(int)neutype];
                  } else {
                    s_CC_data[target][neutype][flavor][e1] = n_sigma_CC_raw_data[e1][1+2*((int)flavor)+(int)neutype];
                    s_NC_data[target][neutype][flavor][e1] = n_sigma_NC_raw_data[e1][1+2*((int)flavor)+(int)neutype];
                  }
                  for (unsigned int e2 = 0; e2 < data_e_size; e2++){
                    if(target == Proton){
                      dsde_CC_data[target][neutype][flavor][e1][e2] = p_dsde_CC_raw_data[e1*data_e_size+e2][2+2*(static_cast<int>(flavor))+static_cast<int>(neutype)];
                      dsde_NC_data[target][neutype][flavor][e1][e2] = p_dsde_NC_raw_data[e1*data_e_size+e2][2+2*(static_cast<int>(flavor))+static_cast<int>(neutype)];
                    } else {
                      dsde_CC_data[target][neutype][flavor][e1][e2] = n_dsde_CC_raw_data[e1*data_e_size+e2][2+2*(static_cast<int>(flavor))+static_cast<int>(neutype)];
                      dsde_NC_data[target][neutype][flavor][e1][e2] = n_dsde_NC_raw_data[e1*data_e_size+e2][2+2*(static_cast<int>(flavor))+static_cast<int>(neutype)];
                    }
                  }
                }
              }
            }
          }
  }
  // declare the object as initialized;
  is_init = true;
}

NeutrinoDISCrossSectionsFromTables::NeutrinoDISCrossSectionsFromTables():
  NeutrinoDISCrossSectionsFromTables(XSECTION_LOCATION "csms.h5"){
}

NeutrinoDISCrossSectionsFromTables::NeutrinoDISCrossSectionsFromTables(std::string path){
  bool is_isoscalar;
  //If a single file, read HDF5
  if(fexists(path)){
    {
      H5File h5file(H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
      readH5Attribute(h5file, "Emin", Emin);
      readH5Attribute(h5file, "Emax", Emax);
      int rank = getArrayH5Rank(h5file, "s_CC");
      if(rank == 4){ // non-isoscalar
        is_isoscalar = false;
        registered_targets_ = {Proton, Neutron};
        readArrayH5(h5file, "s_CC", s_CC_data);
        readArrayH5(h5file, "s_NC", s_NC_data);
        readArrayH5(h5file, "dsDE_CC", dsde_CC_data);
        readArrayH5(h5file, "dsDE_NC", dsde_NC_data);
      }
      if(rank == 3){ // isoscalar
        is_isoscalar = true;
         registered_targets_ = {Isoscalar};

        marray<double,3> s_tmp_data;
        readArrayH5(h5file, "s_CC", s_tmp_data);
        writeIntoLargerMArray<3>(s_tmp_data,s_CC_data);
        readArrayH5(h5file, "s_NC", s_tmp_data);
        writeIntoLargerMArray<3>(s_tmp_data,s_NC_data);

        marray<double,4> dsde_tmp_data;
        readArrayH5(h5file, "dsDE_CC", dsde_CC_data);
        writeIntoLargerMArray<4>(dsde_tmp_data,dsde_CC_data);
        readArrayH5(h5file, "dsDE_NC", dsde_tmp_data);
        writeIntoLargerMArray<4>(dsde_tmp_data,dsde_NC_data);
      }
    }
    //TODO: make sanity checking more user friendly
    int off_set = (is_isoscalar) ? 0 : 1;
    assert(dsde_CC_data.extent(2+off_set)==dsde_CC_data.extent(3+off_set));
    assert(dsde_NC_data.extent(2+off_set)==dsde_NC_data.extent(3+off_set));
    assert(dsde_CC_data.extent(2+off_set)==dsde_NC_data.extent(2+off_set));
    assert(dsde_CC_data.extent(0+off_set)==2);
    assert(dsde_NC_data.extent(0+off_set)==2);
    assert(dsde_CC_data.extent(1+off_set)==3);
    assert(dsde_NC_data.extent(1+off_set)==3);
    
    div=dsde_CC_data.extent(2+off_set);
    unsigned int data_e_size=div;
    auto raw_logE_data_range=logspace(Emin/GeV,Emax/GeV,div);
    logE_data_range.resize(data_e_size);
    std::transform(raw_logE_data_range.begin(),raw_logE_data_range.end(),logE_data_range.begin(),
                   (double(*)(double))log);
    is_init=true;
  }
  else //otherwise, try to read several text files
    ReadText(path);
}

NeutrinoDISCrossSectionsFromTables::~NeutrinoDISCrossSectionsFromTables(){}
    
void NeutrinoDISCrossSectionsFromTables::WriteHDF(std::string path) const{
  H5File h5file(H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
  addH5Attribute(h5file, "Emin", Emin);
  addH5Attribute(h5file, "Emax", Emax);
  writeArrayH5(h5file, "s_CC", s_CC_data, 0);
  writeArrayH5(h5file, "s_NC", s_NC_data, 0);
  writeArrayH5(h5file, "dsDE_CC", dsde_CC_data, 0);
  writeArrayH5(h5file, "dsDE_NC", dsde_NC_data, 0);
}


void NeutrinoDISCrossSectionsFromTables::WriteText(std::string basePath) const{
  std::string filename_sigma_CC = basePath+"sigma_CC.dat";
  std::string filename_sigma_NC = basePath+"sigma_NC.dat";
  std::string filename_dsde_CC = basePath+"dsde_CC.dat";
  std::string filename_dsde_NC = basePath+"dsde_NC.dat";
  
  auto writeTotal=[](const std::string& path, const marray<double,4>& data,
                     const std::vector<double>& logEnergies, double GeV){
    std::ofstream out(path);
    for(size_t ie=0; ie<logEnergies.size(); ie++){
      out << exp(logEnergies[ie]) << ' ';
      for(size_t fl : {0,1,2}){
        for(size_t pt : {0,1})
          for(size_t tg : {0,1,2}){
            out << data[tg][pt][fl][ie] << ' ';
          }
      }
      out << '\n';
    }
  };
  writeTotal(filename_sigma_CC,s_CC_data,logE_data_range,GeV);
  writeTotal(filename_sigma_NC,s_NC_data,logE_data_range,GeV);
  
  auto writeDifferential=[](const std::string& path, const marray<double,5>& data,
                            const std::vector<double>& logEnergies, double GeV){
    std::ofstream out(path);
    for(size_t ie1=0; ie1<logEnergies.size(); ie1++){
      double e1=exp(logEnergies[ie1]);
      for(size_t ie2=0; ie2<logEnergies.size(); ie2++){
        double e2=exp(logEnergies[ie2]);
        out << e1 << ' ' << e2 << ' ';
        for(size_t fl : {0,1,2}){
          for(size_t pt : {0,1})
            for(size_t tg : {0,1,2}){
              out << data[tg][pt][fl][ie1][ie2] << ' ';
            }
        }
        out << '\n';
      }
    }
  };
  writeDifferential(filename_dsde_CC,dsde_CC_data,logE_data_range,GeV);
  writeDifferential(filename_dsde_NC,dsde_NC_data,logE_data_range,GeV);
}
  
GlashowResonanceCrossSection::GlashowResonanceCrossSection(){
  fermi_scale = pow(constants.GF/constants.cm, 2)*constants.electron_mass/constants.pi;
  M_W = 80.385*constants.GeV;
  W_total = 2.085*constants.GeV;
  registered_targets_ = {Electron};
}
  
GlashowResonanceCrossSection::~GlashowResonanceCrossSection(){}
  
double GlashowResonanceCrossSection::TotalCrossSection(double Enu, NeutrinoFlavor flavor,
 NeutrinoType neutype, Current current, Target target) const{
  // only treat the glashow resonance
  if (not (flavor == electron and neutype == antineutrino and current == GR))
    return 0;
  if ( target != Electron )
    return 0.;
  // calculate total cross-section for nuebar + e- -> numubar + mu-, then divide
  // by the W- -> numubar + mu- branching fraction to get total cross-section
  // see e.g. arxiv:1108.3163v2, Eq. 2.1
  const double &m = constants.electron_mass;
  const double &mu = constants.muon_mass;
  return 2*fermi_scale*Enu*pow(1. - (mu*mu - m*m)/(2*m*Enu), 2)
   / (pow(1. - 2*m*Enu/(M_W*M_W), 2) + pow(W_total/M_W, 2)) / B_Muon / 3;
}
  
double GlashowResonanceCrossSection::SingleDifferentialCrossSection(double E1, double E2, NeutrinoFlavor flavor, NeutrinoType neutype, Current current, Target target) const{
  // only treat the glashow resonance
  if (target != Electron)
    return 0.;
  if (not (flavor == electron and neutype == antineutrino and current == GR))
    return 0;
  if (E2 > E1)
    return 0;
  // differential cross section for leptonic final states only
  // NB: assumes that the branching fractions are identical for all 3 families
  double xl = E2/E1;
  return B_Muon*3*TotalCrossSection(E1, flavor, neutype, current, target)*(xl*xl)/E1*constants.GeV;
}
  
double GlashowResonanceCrossSection::DoubleDifferentialCrossSection(double E, double x, double y, NeutrinoFlavor flavor, NeutrinoType neutype, Current current, Target target) const{
  return 0;
}
  
//K. Olive et al. (PDG), Chin. Phys. C38, 090001 (2014)
double GlashowResonanceCrossSection::B_Electron = .1071; //unc. .0016
double GlashowResonanceCrossSection::B_Muon     = .1063; //unc. .0015
double GlashowResonanceCrossSection::B_Tau      = .1138; //unc. .0021
double GlashowResonanceCrossSection::B_Hadronic = .6741; //unc. .0027

} // close namespace
