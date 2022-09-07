#define ProtonNewTreeMaker_run5387_cxx
#include "ProtonNewTreeMaker_run5387.h"

#include <stdexcept>      // std::out_of_range
#include <vector>
#include <iostream>

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <TLegend.h>
#include <TLegendEntry.h>
#include <TMath.h>
#include <TLine.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>
#include <TProfile2D.h>
//#include <iostream>
#include <fstream>
#include <string>
#include "TCanvas.h" 
#include "TGraph.h" 
#include "TGraphSmooth.h" 
#include "TVectorD.h"
#include "TParameter.h"
#include "TGraphErrors.h"
#include <TProfile3D.h>
#include <TProfile2D.h>
#include "TVector3.h"

#include <TMath.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TRandom2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF2.h>
#include <TH1.h>
#include <Math/Functor.h>
#include <TPolyLine3D.h>
#include <Math/Vector3D.h>
#include <Fit/Fitter.h>

//#include "./cali/dedx_function_r5387.h"
#include "./headers/BasicParameters.h"
#include "./headers/BasicAnaFunc.h"
#include "./headers/util.h"
#include "./headers/BetheBloch.h"

#include <cassert>

using namespace std;
using namespace ROOT::Math;


// define the parametric line equation
void line(double t, const double *p, double &x, double &y, double &z) {
	// a parametric line is define from 6 parameters but 4 are independent
	// x0,y0,z0,z1,y1,z1 which are the coordinates of two points on the line
	// can choose z0 = 0 if line not parallel to x-y plane and z1 = 1;
	x = p[0] + p[1]*t;
	y = p[2] + p[3]*t;
	z = t;
}

bool first = true;

// function Object to be minimized
struct SumDistance2 {
	// the TGraph is a data member of the object
	TGraph2D *fGraph;

	SumDistance2(TGraph2D *g) : fGraph(g) {}

	// calculate distance line-point
	double distance2(double x,double y,double z, const double *p) {
		// distance line point is D= | (xp-x0) cross  ux |
		// where ux is direction of line and x0 is a point in the line (like t = 0)
		XYZVector xp(x,y,z);
		XYZVector x0(p[0], p[2], 0. );
		XYZVector x1(p[0] + p[1], p[2] + p[3], 1. );
		XYZVector u = (x1-x0).Unit();
		double d2 = ((xp-x0).Cross(u)).Mag2();
		return d2;
	}

	// implementation of the function to be minimized
	double operator() (const double *par) {
		assert(fGraph != 0);
		double * x = fGraph->GetX();
		double * y = fGraph->GetY();
		double * z = fGraph->GetZ();
		int npoints = fGraph->GetN();
		double sum = 0;
		for (int i  = 0; i < npoints; ++i) {
			double d = distance2(x[i],y[i],z[i],par);
			sum += d;
		}
		if (first) {
			std::cout << "Total Initial distance square = " << sum << std::endl;
		}
		first = false;
		return sum;
	}

};

bool myComparison(const pair<double,int> &a,const pair<double,int> &b) {
	return a.first<b.first;
}

void ProtonNewTreeMaker_run5387::Loop() {
	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntries();
	std::cout<<"nentries: "<<nentries<<std::endl;

	//Energy loss -------------------------------//
	double mean_kecalo_stop=3.78066e+02;
	double err_mean_kecalo_stop=9.81156e-01;
	double sigma_kecalo_stop=5.60675e+01;
	double err_sigma_kecalo_stop=8.12221e-01;

	double mean_kerange_stop=4.00446e+02;
	double err_mean_kerange_stop=1.07727e+00;
	double sigma_kerange_stop=5.12091e+01;
	double err_sigma_kerange_stop=1.06785e+00;

	double mean_kebeam=4.41392e+02;
	double err_mean_kebeam=5.95236e-01;
	double sigma_kebeam=5.15066e+01;
	double err_sigma_kebeam=4.28244e-01;

	double mean_Elossrange_stop=mean_kebeam-mean_kerange_stop;
	//double mean_Elosscalo_stop=mean_kebeam-mean_kecalo_stop;
	//double mean_Elosscalo_stop=(4.85990e+01)/(1.00263e+00); //from fit result
	double mean_Elosscalo_stop=(4.43214e+01)/(9.61985e-01); //from fit result with beam xy

	//fit result ==========================================================//
	//p0           4.85990e+01   8.01697e-01  -2.71531e-08   8.79542e-12
	//p1          -1.00263e+00   1.67201e-02   1.67201e-02   1.66581e-08
	//fit result ==========================================================//

	//-------------------------------------------//

	//Basic configure ------//
	//BetheBloch BB;
	//BB.SetPdgCode(2212);
	//----------------------//

	//Tree Variables -----------------//
	Int_t RUN=-1;
	Int_t subRUN=-1;
	Int_t evt;
	Float_t ntrklen=-1;
	Float_t trklen=-1;
	Float_t PID=-1;
	Float_t B=999;
	Float_t costheta=-999;
	Float_t mediandedx=-999;
	Float_t endpointdedx=-999;
	Float_t calo=-1;
	Float_t avcalo=-1;

	Float_t st_x=-999;
	Float_t st_y=-999;
	Float_t st_z=-999;
	Float_t end_x=-999;
	Float_t end_y=-999;
	Float_t end_z=-999;

	//Float_t pbdt=-999;
	//Int_t nd=0; //num of daughters




	//Tree Structures --------------------------------------------------------------------------------------------------//
	TString str_out=Form("/dune/data2/users/hyliao/protonana/v09_39_01/Classification/protons_data_mva2_run%d.root",run);

	TFile *hfile =new TFile(str_out.Data(),"RECREATE");
	TTree *tree = new TTree("tr","signal");
	tree->Branch("RUN", &RUN, "RUN/I");
	tree->Branch("subRUN", &subRUN, "subRUN/I");
	tree->Branch("evt", &evt, "evt/I");

	tree->Branch("ntrklen", &ntrklen, "ntrklen/F");
	tree->Branch("trklen", &trklen, "trklen/F");
	tree->Branch("B", &B, "B/F");
	tree->Branch("PID", &PID, "PID/F");
	tree->Branch("costheta", &costheta, "costheta/F");
	tree->Branch("mediandedx", &mediandedx, "mediandedx/F");
	tree->Branch("endpointdedx", &endpointdedx, "endpointdedx/F");
	tree->Branch("calo", &calo, "calo/F");
	tree->Branch("avcalo", &avcalo, "avcalo/F");

	tree->Branch("st_x", &st_x, "st_x/F");
	tree->Branch("st_y", &st_x, "st_y/F");
	tree->Branch("st_z", &st_x, "st_z/F");
	tree->Branch("end_x", &st_x, "end_x/F");
	tree->Branch("end_y", &st_x, "end_y/F");
	tree->Branch("end_z", &st_x, "end_z/F");

	//tree->Branch("pbdt", &pbdt, "pbdt/F");
	//tree->Branch("nd", &nd, "nd/I");




	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<nentries;jentry++) { //evt loop
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		if (jentry%10000==0) std::cout<<jentry<<"/"<<nentries<<std::endl;

		//size_t j=0;
		//if (Cut(jentry) < 0) continue;

		//define good beam trigger -----------------------------------------------------//
		bool IsBeamInstValid=false;
		bool IsBeamInst_nMomenta=false;
		bool IsBeamInst_nTracks=false;
		bool IsBeamOK=false;

		if (beam_inst_valid) IsBeamInstValid=true;
		if (beam_inst_nMomenta==1) IsBeamInst_nMomenta=true;
		if (beam_inst_nTracks==1) IsBeamInst_nTracks=true;	
		if (IsBeamInstValid&&IsBeamInst_nMomenta&&IsBeamInst_nTracks) IsBeamOK=true;

		if (!IsBeamOK) continue; //skip evt if beam not okay

		RUN=run;
		subRUN=subrun;
		evt=event;

		//Incoming candidate event ------------------------------------------------------------------------//	
		bool IsProton=false;
		if (!beam_inst_PDG_candidates->empty()&&IsBeamInstValid) {
			//n_can++;
			for (int h=0; h<(int)beam_inst_PDG_candidates->size(); ++h) { 
				//cout<<"beam_inst_PDG_candidates["<<h<<"]:"<<beam_inst_PDG_candidates->at(h)<<endl;
				if (beam_inst_PDG_candidates->at(h)==2212) {
					//n_proton++;
					IsProton=true;
				}
			}
		}
		if (!IsProton) continue; //only protons, not other candidate particles
		//(the ntuple actually only save proton evt for the beam track)

		//Pandora slice cut ------------------------------------------------------------------//
		bool IsPandoraSlice=false; //pandora slice
		bool IsCaloSize=false; //if calo size not empty

		if (isprimarytrack==1&&isprimaryshower==0) { //pandora slice cut
			IsPandoraSlice=true;
		}

		if (!primtrk_hitz->empty()) { 
			IsCaloSize=true;
		}

		//reco pos info & cut --------------------------------------------------------//
		double reco_stx=-99, reco_sty=-99, reco_stz=-99;
		double reco_endx=-99, reco_endy=-99, reco_endz=-99;
		bool IsPos=false;
		if (IsCaloSize) { //calosz
			if (primtrk_hitz->at(0)>primtrk_hitz->at(primtrk_hitz->size()-1)) {
				reco_endx=primtrk_hitx->at(0);
				reco_endy=primtrk_hity->at(0);
				reco_endz=primtrk_hitz->at(0);

				reco_stx=primtrk_hitx->at(primtrk_dedx->size()-1);
				reco_sty=primtrk_hity->at(primtrk_dedx->size()-1);
				reco_stz=primtrk_hitz->at(primtrk_dedx->size()-1);
			}
			else {
				reco_stx=primtrk_hitx->at(0);
				reco_sty=primtrk_hity->at(0);
				reco_stz=primtrk_hitz->at(0);

				reco_endx=primtrk_hitx->at(primtrk_dedx->size()-1);
				reco_endy=primtrk_hity->at(primtrk_dedx->size()-1);
				reco_endz=primtrk_hitz->at(primtrk_dedx->size()-1);
			}

			double beam_dx=(reco_stx-mean_StartX)/sigma_StartX;
			double beam_dy=(reco_sty-mean_StartY)/sigma_StartY;
			double beam_dz=(reco_stz-mean_StartZ)/sigma_StartZ;
			double beam_dxy=sqrt(pow(beam_dx,2)+pow(beam_dy,2));

			if (beam_dx>=dx_min&&beam_dx<=dx_max) { //dx
				if (beam_dy>=dy_min&&beam_dy<=dy_max) { //dy
					if (beam_dz>=dz_min&&beam_dz<=dz_max) { //dz
						if (beam_dxy>=dxy_min&&beam_dxy<=dxy_max) { //dxy
							IsPos=true;
						} //dxy
					} //dz
				} //dy
			} //dx
		} //calosz

		//cosine_theta/cut ---------------------------------------------------------------------------------//
		bool IsCosine=false;
		double cosine_beam_spec_primtrk=-999;
		//cosine_beam_spec_primtrk=cosine_beam_primtrk; //cosine between beam_spec and primary trk direction
		TVector3 dir;
		if (IsCaloSize) {
			//trk direction after SCE corr.
			TVector3 pt0(primtrk_hitx->at(0), primtrk_hity->at(0), primtrk_hitz->at(0));
			TVector3 pt1(primtrk_hitx->at(-1+primtrk_hitx->size()), primtrk_hity->at(-1+primtrk_hity->size()), primtrk_hitz->at(-1+primtrk_hitz->size()));
			dir = pt1 - pt0;
			dir = dir.Unit();

			//beam direction
			//TVector3 beamdir(cos(beam_angleX_mc*TMath::Pi()/180), cos(beam_angleY_mc*TMath::Pi()/180), cos(beam_angleZ_mc*TMath::Pi()/180));
			TVector3 beamdir(beamDirx->at(0), beamDiry->at(0), beamDirz->at(0));
			beamdir = beamdir.Unit();
			cosine_beam_spec_primtrk=dir.Dot(beamdir);
		}
		if (cosine_beam_spec_primtrk<0) { cosine_beam_spec_primtrk=-1.*cosine_beam_spec_primtrk; }
		//if (cosine_beam_spec_primtrk>min_cosine) { IsCosine=true; }
		if (cosine_beam_spec_primtrk>costh_min&&cosine_beam_spec_primtrk<costh_max) { IsCosine=true; }		

		//beam quality cut -----------------------------//
		bool IsBQ=false;
		if (IsCosine&&IsPos) IsBQ=true;
		//if (IsBQ&&IsCaloSize&&IsPandoraSlice) n_bq++;

		//beam info //
		double bx=beamPosx->at(0);
		double by=beamPosy->at(0);
		double bz=beamPosz->at(0);
		double p_beam=beamMomentum->at(0);
		double ke_beam_MeV=1000.*p2ke(p_beam); //unit:MeV
		bool IsBeamMom=false; //apply 3-sigma cut and remove tail events
		if (ke_beam_MeV>=(mean_kebeam-3.*sigma_kebeam)&&ke_beam_MeV<=(mean_kebeam+3.*sigma_kebeam)) IsBeamMom=true;
		if (IsBeamMom==false) continue; //if not within 3-sigma, skip the event

		bool IsBeamXY=false;
		if ((pow(((bx-meanX_data)/(1.5*rmsX_data)),2)+pow(((by-meanY_data)/(1.5*rmsY_data)),2))<=1.) IsBeamXY=true;

		//Get calo info -----------------------------------------------------------------------------------------------//
		double range_reco=-99;
		vector<double> reco_trklen_accum;
		reco_trklen_accum.reserve(primtrk_hitz->size());
		vector<double> EDept;
		double pid=-99;

		vector<double> trkdedx;
		vector<double> trkres;
		double ke_calo_MeV=0;
		double median_dedx=-99;
		double endpoint_dedx=-99;

		//reco tracks
		vector< pair<double,int > > zreco_rawindex; //z, original_index
		//vector< pair<double,double > > zreco_kereco; //z, ke_reco
		vector< pair<double,double > > zreco_zreco; //z, z
		vector< pair<double,double > > zreco_edept; //z, wid_reco
		//vector< pair<double,double > > zreco_widreco; //z, wid_reco
		vector< pair<double,double > > zreco_xreco; //z, wid_reco
		vector< pair<double,double > > zreco_yreco; //z, wid_reco

		vector<double> KE_RECO;
		vector<double> EDEPT_RECO;
		vector<double> WID_RECO;
		vector<double> Z_RECO;
		vector<double> Y_RECO;
		vector<double> X_RECO;
		if (IsCaloSize) { //if calo size not empty
			for (size_t h=0; h<primtrk_dedx->size(); ++h) { //loop over reco hits in a track
				double hitx_reco=primtrk_hitx->at(h);
				double hity_reco=primtrk_hity->at(h);
				double hitz_reco=primtrk_hitz->at(h);
				double resrange_reco=primtrk_resrange->at(h);

				double dqdx=primtrk_dqdx->at(h);
				double dedx=primtrk_dedx->at(h); //prod4-reco 2 has dedx already
				double pitch=primtrk_pitch->at(h);

				int wid_reco=primtrk_wid->at(-1+primtrk_wid->size()-h);
				double pt_reco=primtrk_pt->at(-1+primtrk_wid->size()-h);

				//double cali_dedx=0.;
				//cali_dedx=dedx_function_35ms(dqdx, hitx_reco, hity_reco, hitz_reco);
				EDept.push_back(dedx*pitch);
				ke_calo_MeV+=dedx*pitch;

				if (h==1) range_reco=0;
				if (h>=1) {
					range_reco += sqrt( pow(primtrk_hitx->at(h)-primtrk_hitx->at(h-1), 2)+
							pow(primtrk_hity->at(h)-primtrk_hity->at(h-1), 2)+
							pow(primtrk_hitz->at(h)-primtrk_hitz->at(h-1), 2) );
					reco_trklen_accum[h] = range_reco;
				}

				trkdedx.push_back(dedx);
				trkres.push_back(resrange_reco);

				zreco_rawindex.push_back(make_pair(hitz_reco, h));
				zreco_zreco.push_back(make_pair(hitz_reco, hitz_reco));
				zreco_edept.push_back(make_pair(hitz_reco, dedx*pitch));
				zreco_xreco.push_back(make_pair(hitz_reco, hitx_reco));
				zreco_yreco.push_back(make_pair(hitz_reco, hity_reco));

			} //loop over reco hits in a track
			//range_reco=primtrk_range->at(0); 

			pid=chi2pid(trkdedx,trkres); //pid using stopping proton hypothesis
			median_dedx=TMath::Median(trkdedx.size(), &trkdedx.at(0));
			if (trkdedx.size()>=2) {
				endpoint_dedx=0.5*(trkdedx.at(-1+trkdedx.size())+trkdedx.at(-2+trkdedx.size()));
			}
		} //if calo size not empty	

		//sorting in z --------------------------------------------------------------------------------------//
		sort(zreco_rawindex.begin(),zreco_rawindex.end(),myComparison); //sorting based on the first column
		sort(zreco_zreco.begin(),zreco_zreco.end(),myComparison); //sorting based on the first column
		sort(zreco_edept.begin(),zreco_edept.end(),myComparison); //sorting based on the first column
		sort(zreco_xreco.begin(),zreco_xreco.end(),myComparison); //sorting based on the first column
		sort(zreco_yreco.begin(),zreco_yreco.end(),myComparison); //sorting based on the first column

		for (int si=0; si<(int)zreco_rawindex.size(); ++si){ //save sorting result
			//WID_RECO.push_back(zreco_widreco[si].second);
			//KE_RECO.push_back(zreco_kereco[si].second);
			EDEPT_RECO.push_back(zreco_edept[si].second);
			Z_RECO.push_back(zreco_zreco[si].first);
			Y_RECO.push_back(zreco_yreco[si].second);
			X_RECO.push_back(zreco_xreco[si].second);
		} //save sorting result

		//cout<<"Z_RECO.size="<<Z_RECO.size()<<endl;
		//b2 calculation using 3D fit -------------------------------------------------------------------------------------------------------------------------//
		TGraph2D *gr=new TGraph2D();
		if (Z_RECO.size()>2) { //at leat 3 hits
			int fit_zmin=0; //start hit for fitting
			int fit_zmax=24; //end fit for fitting
			if ((int)Z_RECO.size()<=fit_zmax) fit_zmax=-1+(int)Z_RECO.size();

			for (int N=fit_zmin; N<fit_zmax; N++) {
				gr->SetPoint(N,X_RECO.at(N), Y_RECO.at(N), Z_RECO.at(N));
			}

			//Initialization of parameters for 3d fit
			int N=(int)Z_RECO.size();
			double ini_p1=(X_RECO.at(N-1)-X_RECO.at(0))/(Z_RECO.at(N-1)-Z_RECO.at(0));
			double ini_p0=X_RECO.at(0)-ini_p1*Z_RECO.at(0);
			double ini_p3=Y_RECO.at(N-1)-Y_RECO.at(0);
			double ini_p2=Y_RECO.at(0)-ini_p3*Z_RECO.at(0);

			ROOT::Fit::Fitter  fitter;
			// make the functor objet
			SumDistance2 sdist(gr);
			ROOT::Math::Functor fcn(sdist,4);

			// set the function and the initial parameter values
			double pStart[4]={ini_p0, ini_p1, ini_p2, ini_p3};   
			fitter.SetFCN(fcn,pStart);

			// set step sizes different than default ones (0.3 times parameter values)
			for (int ik = 0; ik < 4; ++ik) fitter.Config().ParSettings(ik).SetStepSize(0.01);

			bool ok = fitter.FitFCN();
			if (!ok) {
				//Error("line3Dfit","Line3D Fit failed");
				//return 1;
			}

			const ROOT::Fit::FitResult & result = fitter.Result();
			//std::cout << "Total final distance square " << result.MinFcnValue() << std::endl;
			//result.Print(std::cout);

			// get fit parameters
			const double * parFit = result.GetParams();

			//cross-product to calculate b2
			//         A (end of track)
			//         *\
			//         | \
			//         |  \
			//         |   *B (start of track)
			// *------------
			// C (selected end of line)

			//start-point in the line
			double zl_st=Z_RECO.at(0);
			double yl_st=result.Parameter(2)+result.Parameter(3)*zl_st;
			double xl_st=result.Parameter(0)+result.Parameter(1)*zl_st;

			//selected end-point in the line
			double zl_end=Z_RECO.at(N-1);
			double yl_end=result.Parameter(2)+result.Parameter(3)*zl_end;
			double xl_end=result.Parameter(0)+result.Parameter(1)*zl_end;

			//b2 calculation
			TVector3 BC;
			BC.SetXYZ(xl_end-xl_st, yl_end-yl_st, zl_end-zl_st);

			TVector3 BA;
			BA.SetXYZ(X_RECO.at(N-1)-xl_st, Y_RECO.at(N-1)-yl_st, Z_RECO.at(N-1)-zl_st);

			B=(BA.Cross(BC)).Mag()/BC.Mag();
			//std::cout<<"Minimum distance (b2):"<<B<<std::endl;
		} //at leat 3 hits

		//clean memory
		delete gr;

		//b2 calculation using 3D fit -------------------------------------------------------------------------------------------------------------------------//




		double csda_val=csda_range_vs_mom_sm->Eval(p_beam); //expected csda value if proton stops; unit: cm
		double ke_range=ke_vs_csda_range_sm->Eval(range_reco); //unit:GeV
		double ke_range_MeV=1000.*ke_range; //unit:MeV
		double p_range_MeV=1000.*ke2p(ke_range_MeV/1000.); //MeV/c
		double norm_trklen=range_reco/csda_val; //trklen/csda_val
		//cout<<"csda_val:"<<csda_val<<" range_reco:"<<range_reco<<" norm_trklen:"<<norm_trklen<<" min_norm_trklen_csda:"<<min_norm_trklen_csda<<endl;

		bool IsRecoStop=false;
		bool IsRecoInEL=false;
		bool IsRecoEL=false;
		if (norm_trklen>min_norm_trklen_csda) IsRecoStop=true; //stopping proton cut
		//if (norm_trklen<=min_norm_trklen_csda) IsRecoInel=true; //inelastic proton cut
		//if (IsRecoInel&&IsBQ&&IsCaloSize&&IsPandoraSlice) n_recoinel++;

		if ((range_reco/csda_val)<min_norm_trklen_csda) { //inel region
			if (pid>pid_1) IsRecoInEL=true;
			if (pid<=pid_1) IsRecoEL=true;
		} //inel region
		if ((range_reco/csda_val)>=min_norm_trklen_csda&&(range_reco/csda_val)<max_norm_trklen_csda) { //stopping p region
			if (pid>pid_2) IsRecoInEL=true;
			if (pid<=pid_2) IsRecoEL=true;
		} //stopping p region

		//XY Cut ---------------------------------------------------------------------------------------------------------------//
		bool IsXY=false;
		//start(x,y,z) without SCE corr. 
		double xst_nosce=0;
		double yst_nosce=0;
		double zst_nosce=0;
		if (IsCaloSize&&IsPandoraSlice) { //if calo size not empty
			if ((primtrk_startz->at(-1+primtrk_startz->size()))>(primtrk_startz->at(0))) { //check if Pandora flip the sign
				xst_nosce=primtrk_startx->at(0);
				yst_nosce=primtrk_starty->at(0);
				zst_nosce=primtrk_startz->at(0);
			} //check if Pandora flip the sign
			else {
				xst_nosce=primtrk_startx->at(-1+primtrk_startx->size());
				yst_nosce=primtrk_starty->at(-1+primtrk_starty->size());
				zst_nosce=primtrk_startz->at(-1+primtrk_startz->size());
			}
		} //if calo size not empty
		if ((pow(((xst_nosce-mean_x)/dev_x),2)+pow(((yst_nosce-mean_y)/dev_y),2))<=1.) IsXY=true;
		//-----------------------------------------------------------------------//

		bool IsIntersection=false;
		if (timeintersection->size()) IsIntersection=true;

		//calo info
		//double ke_calo_MeV=0;
		double xst_sce=0;
		double yst_sce=0;
		double zst_sce=0;
		double xend_sce=0;
		double yend_sce=0;
		double zend_sce=0;
		if (IsBQ&&IsCaloSize&&IsPandoraSlice) { //if calo size not empty
			//h2d_xy_noSCE->Fill(xst_nosce, yst_nosce);
			//h1d_zst_noSCE->Fill(zst_nosce);

			//Fill1DHist(h1d_KEbb_reco, KEbb_reco2);
			//Fill1DHist(h1d_KEcalo, keff_reco4-ke_calo_MeV);

			//start(x,y,z) after SCE corr. ----------------------------------------------------------------------------//
			if ((primtrk_hitz->at(-1+primtrk_hitz->size()))>(primtrk_hitz->at(0))) { //check if Pandora flip the sign
				xst_sce=primtrk_hitx->at(0);
				yst_sce=primtrk_hity->at(0);
				zst_sce=primtrk_hitz->at(0);

				xend_sce=primtrk_hitx->at(-1+primtrk_hitx->size());
				yend_sce=primtrk_hity->at(-1+primtrk_hity->size());
				zend_sce=primtrk_hitz->at(-1+primtrk_hitz->size());
			} //check if Pandora flip the sign
			else {
				xst_sce=primtrk_hitx->at(-1+primtrk_hitx->size());
				yst_sce=primtrk_hity->at(-1+primtrk_hity->size());
				zst_sce=primtrk_hitz->at(-1+primtrk_hitz->size());

				xend_sce=primtrk_hitx->at(0);
				yend_sce=primtrk_hity->at(0);
				zend_sce=primtrk_hitz->at(0);
			}
			//h2d_xy_SCE->Fill(xst_sce,yst_sce);
			//h1d_zst_SCE->Fill(zst_sce);


			//ntrklen_chi2pid_BQ->Fill(norm_trklen, pid);
			//Fill1DHist(h1d_chi2pid_BQ, pid);
			//Fill1DHist(h1d_ntrklen_BQ, norm_trklen);

			//Fill1DHist(h1d_mediandedx_BQ, median_dedx);
			//Fill1DHist(h1d_dEdL_BQ, ke_calo_MeV/range_reco);			

			//if (IsRecoStop) { //reco_stop
			//chi2pid_recostop->Fill(pid);
			//} //reco_stop	
			//if (IsRecoInel) { //reco_inel
			//chi2pid_recoinel->Fill(pid);
			//} //reco_inel
		} //if calo size not empty

		//KEff, KEend -----------------------------------------------------------------------------------------//
		double KEcsda=1000.*ke_vs_csda_range_sm->Eval(range_reco);
		double KEcsda_inel=ke_beam_MeV-mean_Elossrange_stop;

		//double KEbb_reco_constErange=-1; KEbb_reco_constErange=BB.KEAtLength(KEcsda_inel, range_reco);
		//double KEbb_reco_Edept_range=-1; KEbb_reco_Edept_range=BB.KEAtLength(KEcsda, range_reco);

		double KEcalo_reco_constEcalo=-1; KEcalo_reco_constEcalo=ke_beam_MeV-mean_Elosscalo_stop-ke_calo_MeV;
		double KEff=ke_beam_MeV-mean_Elosscalo_stop;


		//hypothetical length -------------------------------------------------------------------------------------//
		double fitted_length=-1;
		std::reverse(trkdedx.begin(),trkdedx.end());  
		std::reverse(trkres.begin(),trkres.end());  
		//double tmp_fitted_length=BB.Fit_dEdx_Residual_Length(trkdedx, trkres, 2212, false);
		//if (tmp_fitted_length>0) fitted_length=tmp_fitted_length;
		//double fitted_KE=-50; 
		//if (fitted_length>0) fitted_KE=BB.KEFromRangeSpline(fitted_length);

		//if (IsBQ&&IsCaloSize&&IsPandoraSlice) {
		//if (IsIntersection&&IsBQ&&IsCaloSize&&IsPandoraSlice) {
		if (IsBeamXY&&IsBQ&&IsCaloSize&&IsPandoraSlice) { //basic cuts
			ntrklen=norm_trklen;
			trklen=range_reco;
			PID=pid;
			mediandedx=median_dedx;
			endpointdedx=endpoint_dedx;
			calo=ke_calo_MeV;
			avcalo=ke_calo_MeV/range_reco;

			st_x=reco_stx;
			st_y=reco_sty;
			st_z=reco_stz;

			end_x=reco_endx;
			end_y=reco_endy;
			end_z=reco_endz;

			costheta=cosine_beam_spec_primtrk;

			tree->Fill();			
		} //basic cuts


	} //evt loop


	//save results -------//
	tree->Write();


	}
