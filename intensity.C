// Based on May's thesis, pages 153 and 154.



const double pi=TMath::Pi();
const double small_gamma=30.; // Hz/uT
const double small_gamma_bar=2*pi*small_gamma; // rad/s/uT
const double b0=1; // uT
const double omega_r=small_gamma_bar*b0; // rad/s
const double t_precess=100; // s
const double cap_gamma=pi/t_precess; // rad/s


double fit_up(double *x, double *par) {
  return par[0]+par[1]*pow(sin(pi/2/cap_gamma*(x[0]-par[2])),2);
}

double fit_down(double *x, double *par) {
  return par[0]+par[1]*pow(cos(pi/2/cap_gamma*(x[0]-par[2])),2);
}

void intensity() {

  double alpha=0.8; // dimensionless
  double n_ucn=20000; // total number of neutrons stored

  double n_up_max=n_ucn*alpha;
  double n_up_min=n_ucn*(1-alpha);
  double n_down_max=n_ucn*alpha;
  double n_down_min=n_ucn*(1-alpha);

  double cap_s=0.3; // see the thesis

  const int m=4; // number of measurements to conduct
  double omega[m];
  double n_up[m],n_down[m];
  double n_up_random[m],n_down_random[m];
  double dn_up_random[m],dn_down_random[m];
  double domega[m];

  TRandom *rd = new TRandom();

  for (int i=0;i<m;i++) {
    if (i==0) omega[i]=omega_r+cap_gamma/2*(1+cap_s);
    else if (i==1) omega[i]=omega_r-cap_gamma/2*(1+cap_s);
    else if (i==2) omega[i]=omega_r+cap_gamma/2*(1-cap_s);
    else omega[i]=omega_r-cap_gamma/2*(1-cap_s);
    cout << "omega[i] = " << omega[i] << endl;

    // Calculate expected number UCN, spin up and down
    // May eqs. 5.9 and 5.10
    n_up[i]=n_up_min
      +(n_up_max-n_up_min)*pow(sin(pi/2/cap_gamma*(omega[i]-omega_r)),2);
    n_down[i]=n_down_min
      +(n_down_max-n_down_min)*pow(cos(pi/2/cap_gamma*(omega[i]-omega_r)),2);
    cout << "n_up[i] = " << n_up[i] << endl;
    cout << "n_down[i] = " << n_down[i] << endl;

    // Poisson distribute based on these means
    n_up_random[i]=rd->PoissonD(n_up[i]);
    n_down_random[i]=rd->PoissonD(n_down[i]);

    cout << "n_up_random[i] = " << n_up_random[i] << endl;
    cout << "n_down_random[i] = " << n_down_random[i] << endl;

    // fill errors

    dn_up_random[i]=sqrt(n_up_random[i]);
    dn_down_random[i]=sqrt(n_down_random[i]);
    domega[i]=0;

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

  // Get the most interesting fit parameters
  double omega_up_fit,domega_up_fit;
  omega_up_fit=fitfunc_up->GetParameter(2);
  domega_up_fit=fitfunc_up->GetParError(2);
  cout << omega_up_fit << " " << domega_up_fit << endl;
  double omega_down_fit,domega_down_fit;
  omega_down_fit=fitfunc_down->GetParameter(2);
  domega_down_fit=fitfunc_down->GetParError(2);
  cout << omega_down_fit << " " << domega_down_fit << endl;

  // Compare to equation 5.12 of May's thesis
  cout << cap_gamma/pi/alpha/sqrt(n_ucn) << endl;

  c1->Update();
}
