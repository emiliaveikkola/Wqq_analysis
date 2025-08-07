#include <iostream>

double weight(int flavor, int tag) {

  double h3all = 1;
  double h3gencsall = 0.289841;
  double h3genudall = 0.377007;
  double h3genxall = 0.333153;
  double h3tagcsall = 0.117846;
  double h3tagudall = 0.180028;
  double h3tagxall = 0.702127;
  double h3gencstagcs = 0.0820791;
  double h3gencstagud = 0.0243256;
  double h3gencstagx = 0.183436;
  double h3genudtagcs = 0.00609975;
  double h3genudtagud = 0.091527;
  double h3genudtagx = 0.27938;
  double h3genxtagcs = 0.0296669;
  double h3genxtagud = 0.0641751;
  double h3genxtagx = 0.23931;
  double h3all_data = 1;
  double h3tagcsall_data = 0.0952018;
  double h3tagudall_data = 0.148134;
  double h3tagxall_data = 0.756665;
    /*
    double h3all = 0.117564;
  double h3gencsall = 0.0340747;
  double h3genudall = 0.0443223;
  double h3genxall = 0.0391667;
  double h3tagcsall = 0.0138544;
  double h3tagudall = 0.0211647;
  double h3tagxall = 0.0825446;
  double h3gencstagcs = 0.00964952;
  double h3gencstagud = 0.00285981;
  double h3gencstagx = 0.0215654;
  double h3genudtagcs = 0.00071711;
  double h3genudtagud = 0.0107603;
  double h3genudtagx = 0.032845;
  double h3genxtagcs = 0.00348776;
  double h3genxtagud = 0.00754467;
  double h3genxtagx = 0.0281342;
  double h3all_data = 1;
  double h3tagcsall_data = 0.0952018;
  double h3tagudall_data = 0.148134;
  double h3tagxall_data = 0.756665;
    */

  // Hard-coded efficiencies from MC simulations
  const double epsilon_c_cs = h3gencstagcs/h3gencsall;//0.8;  // Efficiency for cs pair being c-tagged
  const double epsilon_u_cs = h3gencstagud/h3gencsall;//;0.1;  // Efficiency for cs pair being u-tagged
  const double epsilon_c_ud = h3genudtagcs/h3genudall;//0.05; // Efficiency for ud pair being c-tagged
  const double epsilon_u_ud = h3genudtagud/h3genudall;//0.7;  // Efficiency for ud pair being u-tagged
  const double epsilon_c_xx = h3genxtagcs/h3genxall;//0.05; // Efficiency for xx pair being c-tagged
  const double epsilon_u_xx = h3genxtagud/h3genxall;//0.7;  // Efficiency for xx pair being u-tagged
  
  // Event counts from MC and data
  const double N_cs = h3gencsall/h3all;//1000.0;       // Number of cs events
  const double N_ud = h3genudall/h3all;//h1000.0;       // Number of ud events
  const double N_xx = h3genxall/h3all;//h1000.0;       // Number of ud events
  const double N_tag_cs = h3tagcsall_data/h3all_data;///750.0;    // Number of events tagged as cs in data
  const double N_tag_ud = h3tagudall_data/h3all_data;//650.0;    // Number of events tagged as ud in data

  // Assumed scale factors
  const double k_c_ud = 0.8;  // Scale factor for ud pair being c-tagged
  const double k_u_cs = 0.8;  // Scale factor for cs pair being u-tagged
  const double k_c_cs_init = 0.80858; // Initial value for k_c_xx calculation
  const double k_u_ud_init = 0.827996; // Initial value for k_u_xx calculation
  //const double k_c_xx = 1.0;  // Scale factor for xx pair being c-tagged
  //const double k_u_xx = 1.0;  // Scale factor for xx pair being u-tagged

  // Estimate fraction of cs and ud like events in category x
  // epsilon_c_xx = epsilon_c_cs*f_cs_c + epsilon_c_ud*(1-f_cs_c)
  double f_cs_c = (epsilon_c_xx-epsilon_c_ud) / (epsilon_c_cs-epsilon_c_ud);
  // epsilon_u_xx = epsilon_u_cs*f_cs_u + epsilon_u_ud*(1-f_cs_u)
  double f_cs_u = (epsilon_u_xx-epsilon_u_ud) / (epsilon_u_cs-epsilon_u_ud);

  // Weighted efficiency for x category
  const double k_c_xx = (epsilon_c_cs*f_cs_c*k_c_cs_init + epsilon_c_ud*(1-f_cs_c)*k_c_ud) / epsilon_c_xx;
  const double k_u_xx = (epsilon_u_cs*f_cs_u*k_u_cs + epsilon_u_ud*(1-f_cs_u)*k_u_ud_init) / epsilon_u_xx;
  
  // Calculate scale factors based on provided formulas
  double k_c_cs = (N_tag_cs - N_ud * epsilon_c_ud * k_c_ud - N_xx * epsilon_c_xx * k_c_xx) / (N_cs * epsilon_c_cs);
  double k_u_ud = (N_tag_ud - N_cs * epsilon_u_cs * k_u_cs - N_xx * epsilon_u_xx * k_u_xx) / (N_ud * epsilon_u_ud);
  
  // Initialize the weight
  double weight = 1.0;

  // Determine the weight based on flavor and tag
  if (flavor == 4) { // cs pair
    if (tag == 4) { // c-tag
      weight = k_c_cs;
    } else if (tag == 1) { // ud-tag
      weight = k_u_cs;
    } else if (tag == 0) { // xx-tag (untagged)
      double numerator = 1.0 - epsilon_c_cs * k_c_cs - epsilon_u_cs * k_u_cs;
      double denominator = 1.0 - epsilon_c_cs - epsilon_u_cs;
      weight = numerator / denominator;
    } else {
      // Invalid tag category
      weight = 1.0;
    }
  } else if (flavor == 1) { // ud pair
    if (tag == 4) { // c-tag
      weight = k_c_ud;
    } else if (tag == 1) { // ud-tag
      weight = k_u_ud;
    } else if (tag == 0) { // xx-tag (untagged)
      double numerator = 1.0 - epsilon_c_ud * k_c_ud - epsilon_u_ud * k_u_ud;
      double denominator = 1.0 - epsilon_c_ud - epsilon_u_ud;
      weight = numerator / denominator;
    } else {
      // Invalid tag category
      weight = 1.0;
    }
  } else if (flavor == 0) { // xx pair
    if (tag == 4) { // c-tag
      weight = k_c_xx;
    } else if (tag == 1) { // ud-tag
      weight = k_u_xx;
    } else if (tag == 0) { // xx-tag (untagged)
	double numerator = 1.0 - epsilon_c_xx * k_c_xx - epsilon_u_xx * k_u_xx;
	double denominator = 1.0 - epsilon_c_xx - epsilon_u_xx;
	weight = numerator / denominator;
    } else {
      // Invalid tag category
      weight = 1.0;
    }
  } else {
    // Invalid flavor category
    weight = 1.0;
  }
  
  if (bool debug = (flavor==0 && tag==1)) {
    double N_tag_cs_MC = N_cs*epsilon_c_cs + N_ud*epsilon_c_ud + N_xx*epsilon_c_xx;
    double N_tag_ud_MC = N_cs*epsilon_u_cs + N_ud*epsilon_u_ud + N_xx*epsilon_u_xx;
    cout << weight<<"\nflavor="<<flavor<< ", tag="<<tag<<endl;
    cout << "epsilon_c_cs = "<<epsilon_c_cs << endl;
    cout << "epsilon_c_xx = "<<epsilon_c_xx << " (fc="<<f_cs_c<<")"<< endl;
    cout << "epsilon_c_ud = "<<epsilon_c_ud << endl;
    cout << "epsilon_u_cs = "<<epsilon_u_cs << endl;
    cout << "epsilon_u_xx = "<<epsilon_u_xx << " (fc="<<f_cs_u<<")"<< endl;
    cout << "epsilon_u_ud = "<<epsilon_u_ud << endl;
    cout << "N_cs = "<<N_cs << endl;
    cout << "N_ud = "<<N_ud << endl;
    cout << "N_xx = "<<N_xx << endl;
    cout << "N_tag_cs_MC = "<<N_tag_cs_MC << endl;
    cout << "N_tag_cs_DT = "<<N_tag_cs <<" ("<<N_tag_cs/N_tag_cs_MC<<")"<<endl;
    cout << "N_tag_ud_MC = "<<N_tag_ud_MC << endl;
    cout << "N_tag_ud_DT = "<<N_tag_ud <<" ("<<N_tag_ud/N_tag_ud_MC<<")"<<endl;
    cout << "k_c_ud = "<<k_c_ud << endl;
    cout << "k_c_xx = "<<k_c_xx << endl;
    cout << "k_c_cs = "<<k_c_cs << "** (" << k_c_cs_init << ")" << endl;
    cout << "k_u_cs =  "<<k_u_cs << endl;
    cout << "k_u_xx = "<<k_u_xx << endl;
    cout << "k_u_ud = "<<k_u_ud << "** (" << k_u_ud_init << ")" << endl;
  }
  
  return weight;
}

// Example usage
void weight() {
    // Test cases for each category
    std::cout << "Weight for cs pair with c-tag: " << weight(4, 4) << std::endl;
    std::cout << "Weight for cs pair with x-tag: " << weight(4, 0) << std::endl;
    std::cout << "Weight for cs pair with u-tag: " << weight(4, 1) << std::endl;
    std::cout << "Weight for ud pair with c-tag: " << weight(1, 4) << std::endl;
    std::cout << "Weight for ud pair with x-tag: " << weight(1, 0) << std::endl;
    std::cout << "Weight for ud pair with u-tag: " << weight(1, 1) << std::endl;
    std::cout << "Weight for xx pair with c-tag: " << weight(0, 4) << std::endl;
    std::cout << "Weight for xx pair with x-tag: " << weight(0, 0) << std::endl;
    std::cout << "Weight for xx pair with u-tag: " << weight(0, 1) << std::endl;
}
