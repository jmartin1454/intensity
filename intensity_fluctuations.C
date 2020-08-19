const double pi=TMath::Pi();
const double small_gamma=30.; // Hz/uT
const double small_gamma_bar=2*pi*small_gamma; // rad/s/uT
const double b0=1; // uT
const double omega_r=small_gamma_bar*b0; // Hz
const double t_precess=100; // s
const double cap_gamma=pi/t_precess; // Hz


double fit_up(double *x, double *par) {
  return par[0]+par[1]*pow(sin(pi/2/cap_gamma*(x[0]-par[2])),2);
}

double fit_down(double *x, double *par) {
  return par[0]+par[1]*pow(cos(pi/2/cap_gamma*(x[0]-par[2])),2);
}

double fluctuate(TRandom *rd) {
  // returns a random number, Gaussian distributed about 1
  double mean=1;
  //  double sigma=0.04;
  double sigma=0.;
  return rd->Gaus(mean,sigma);
}

void intensity_fluctuations() {

  // Based on May's thesis, pages 153 and 154.

  double alpha=0.5; // dimensionless
  double n_ucn=20000; // total number of neutrons stored on average

  double n_ucn_fluct; // total number stored, with fluctuation
  double n_up_max;
  double n_up_min;
  double n_down_max;
  double n_down_min;

  double cap_s=0.3; // see the thesis

  const int m=4; // number of measurements to conduct
  double omega[m];
  double n_up[m],n_down[m];
  double n_up_random[m],n_down_random[m];
  double dn_up_random[m],dn_down_random[m];
  double domega[m];

  double a[m],da[m];

  TRandom *rd = new TRandom();

  for (int i=0;i<m;i++) {
    // select applied RF frequency for this fill
    if (i==0) omega[i]=omega_r+cap_gamma/2*(1+cap_s);
    else if (i==1) omega[i]=omega_r-cap_gamma/2*(1+cap_s);
    else if (i==2) omega[i]=omega_r+cap_gamma/2*(1-cap_s);
    else omega[i]=omega_r-cap_gamma/2*(1-cap_s);
    cout << "omega[i] = " << omega[i] << endl;

    // generate random intensity fluctuation
    n_ucn_fluct=n_ucn*fluctuate(rd);

    cout << "average Total number of UCN for this fill "
	 << n_ucn_fluct << endl;


    // Calculate expected number UCN, spin up and down
    // May eqs. 5.9 and 5.10
    n_up_max=n_ucn_fluct*(1+alpha)/2;
    n_up_min=n_ucn_fluct*(1-alpha)/2;
    n_down_max=n_ucn_fluct*(1+alpha)/2;
    n_down_min=n_ucn_fluct*(1-alpha)/2;

    n_up[i]=n_up_min
      +(n_up_max-n_up_min)*pow(sin(pi/2/cap_gamma*(omega[i]-omega_r)),2);
    n_down[i]=n_down_min
      +(n_down_max-n_down_min)*pow(cos(pi/2/cap_gamma*(omega[i]-omega_r)),2);
    cout << "n_up[i] = " << n_up[i] << endl;
    cout << "n_down[i] = " << n_down[i] << endl;

    // Poisson distribute based on these as means
    n_up_random[i]=rd->PoissonD(n_up[i]);
    n_down_random[i]=rd->PoissonD(n_down[i]);

    cout << "n_up_random[i] = " << n_up_random[i] << endl;
    cout << "n_down_random[i] = " << n_down_random[i] << endl;

    // fill error arrays
    dn_up_random[i]=sqrt(n_up_random[i]);
    dn_down_random[i]=sqrt(n_down_random[i]);
    domega[i]=0;

    // fill "asymmetry" arrays
    a[i]=(n_up_random[i]-n_down_random[i])/(n_up_random[i]+n_down_random[i]);
    da[i]=sqrt((1-pow(a[i],2))/(n_up_random[i]+n_down_random[i]));
  }

  // Now graph what we got
  TCanvas *c1 = new TCanvas("c1","A Simple Graph with error bars",200,10,700,500);

  TGraphErrors *gr_up =
    new TGraphErrors(m,omega,n_up_random,domega,dn_up_random);
  TGraphErrors *gr_down =
    new TGraphErrors(m,omega,n_down_random,domega,dn_down_random);

  gr_up->SetMarkerStyle(21);
  gr_up->Draw("AP");
  gr_down->SetMarkerStyle(22);
  gr_down->Draw("P");

  // Do the fits
  TF1 *fitfunc_up =
    new TF1("fitfunc_up",fit_up,omega_r-cap_gamma,omega_r+cap_gamma,3);
  fitfunc_up->SetParameters(n_up_min,n_up_max-n_up_min,omega_r);
  gr_up->Fit("fitfunc_up","e");

  TF1 *fitfunc_down =
    new TF1("fitfunc_down",fit_down,omega_r-cap_gamma,omega_r+cap_gamma,3);
  fitfunc_down->SetParameters(n_down_min,n_down_max-n_down_min,omega_r);
  gr_down->Fit("fitfunc_down","e");

  double omega_up=fitfunc_up->GetParameter(2);
  double domega_up=fitfunc_up->GetParError(2);
  double omega_down=fitfunc_down->GetParameter(2);
  double domega_down=fitfunc_down->GetParError(2);

  cout << "May says " << cap_gamma/pi/alpha/sqrt(n_ucn*4) << endl; // May equation 5.11

  // average the measurements together as suggested by May after eq. 5.10
  double denominator=1/pow(domega_up,2)+1/pow(domega_down,2);
  double numerator=omega_up/pow(domega_up,2)+omega_down/pow(domega_down,2);
  double omega_meas=numerator/denominator;
  double domega_meas=sqrt(1/denominator);

  printf("Final answer:  omega = %10.6f +/- %10.6f rad/s\n",
	 omega_meas,domega_meas);

  printf("Correct answer:  omega = %10.6f rad/s\n",omega_r);

  printf("Measured minus correct = %10.6f rad/s\n",omega_meas-omega_r);

  c1->Update();

  TCanvas *c2 = new TCanvas("c2","A Simple Graph with error bars",200,10,700,500);

  TGraphErrors *gr_a = new TGraphErrors(m,omega,a,domega,da);
  gr_a->SetMarkerStyle(23);
  gr_a->Draw("AP");

  c2->Update();

}

