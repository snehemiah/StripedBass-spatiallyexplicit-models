#ifdef DEBUG
  #ifndef __SUNPRO_C
    #include <cfenv>
    #include <cstdlib>
  #endif
#endif
#include <admodel.h>
#include <contrib.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <scaa-stripedbass.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  fdyear.allocate("fdyear");
  ldyear.allocate("ldyear");
  fmyear.allocate("fmyear");
  lmyear.allocate("lmyear");
  fage.allocate("fage");
  lage.allocate("lage");
  nyrs.allocate("nyrs");
  fleet.allocate("fleet");
  region.allocate("region");
  tstep.allocate("tstep");
  stock.allocate("stock");
  nsurv.allocate("nsurv");
  age1surv.allocate("age1surv");
  yoysurv_coast.allocate("yoysurv_coast");
  yoysurv_bay.allocate("yoysurv_bay");
  tblock.allocate("tblock");
  ftyear.allocate(1,tblock,"ftyear");
  lbound.allocate("lbound");
  ubound.allocate("ubound");
  obs_C.allocate(1,region,1,tstep,fdyear,ldyear,"obs_C");
  obs_Cp.allocate(1,region,1,tstep,fdyear,ldyear,fage,lage,"obs_Cp");
  obs_CV.allocate(1,region,1,tstep,fdyear,ldyear,"obs_CV");
  sfage_a.allocate("sfage_a");
  sfage_b.allocate("sfage_b");
  sfage_c.allocate("sfage_c");
  slage_a.allocate("slage_a");
  slage_b.allocate("slage_b");
  slage_c.allocate("slage_c");
  obs_I_md.allocate(fdyear,ldyear,"obs_I_md");
  obs_I_ny.allocate(fdyear,ldyear,"obs_I_ny");
  obs_I_nj.allocate(fdyear,ldyear,"obs_I_nj");
  obs_I_des.allocate(fdyear,ldyear,"obs_I_des");
  obs_I_de30.allocate(fdyear,ldyear,"obs_I_de30");
  obs_I_cm.allocate(1,tstep,fdyear,ldyear,"obs_I_cm");
  obs_I_ct.allocate(1,tstep,fdyear,ldyear,"obs_I_ct");
  obs_Ip_md.allocate(fdyear,ldyear,sfage_b,slage_b,"obs_Ip_md");
  obs_Ip_ny.allocate(fdyear,ldyear,sfage_c,slage_c,"obs_Ip_ny");
  obs_Ip_nj.allocate(fdyear,ldyear,sfage_b,slage_b,"obs_Ip_nj");
  obs_Ip_des.allocate(fdyear,ldyear,sfage_c,slage_c,"obs_Ip_des");
  obs_Ip_de30.allocate(fdyear,ldyear,sfage_a,slage_a,"obs_Ip_de30");
  obs_Ip_cm.allocate(1,tstep,fdyear,ldyear,sfage_a,slage_a,"obs_Ip_cm");
  obs_Ip_ct.allocate(1,tstep,fdyear,ldyear,sfage_a,slage_a,"obs_Ip_ct");
  obs_I_CV_cm.allocate(1,tstep,fdyear,ldyear,"obs_I_CV_cm");
  obs_I_CV_de30.allocate(fdyear,ldyear,"obs_I_CV_de30");
  obs_I_CV_ct.allocate(1,tstep,fdyear,ldyear,"obs_I_CV_ct");
  obs_I_CV_md.allocate(fdyear,ldyear,"obs_I_CV_md");
  obs_I_CV_nj.allocate(fdyear,ldyear,"obs_I_CV_nj");
  obs_I_CV_des.allocate(fdyear,ldyear,"obs_I_CV_des");
  obs_I_CV_ny.allocate(fdyear,ldyear,"obs_I_CV_ny");
  obs_I_age1n.allocate(1,age1surv,fdyear,ldyear,"obs_I_age1n");
  obs_I_age1m.allocate(fdyear,ldyear,"obs_I_age1m");
  obs_I_age1n_CV.allocate(1,age1surv,fdyear,ldyear,"obs_I_age1n_CV");
  obs_I_age1m_CV.allocate(fdyear,ldyear,"obs_I_age1m_CV");
  obs_I_yoy_coast.allocate(1,yoysurv_coast,fdyear,ldyear,"obs_I_yoy_coast");
  obs_I_yoy_bay.allocate(1,yoysurv_bay,fdyear,ldyear,"obs_I_yoy_bay");
  obs_I_yoy_CV_coast.allocate(1,yoysurv_coast,fdyear,ldyear,"obs_I_yoy_CV_coast");
  obs_I_yoy_CV_bay.allocate(1,yoysurv_bay,fdyear,ldyear,"obs_I_yoy_CV_bay");
  w_age.allocate(fdyear,ldyear,fage,lage,"w_age");
  ssbw_age.allocate(fdyear,ldyear,fage,lage,"ssbw_age");
  rw_age.allocate(fdyear,ldyear,fage,lage,"rw_age");
  m_age.allocate(fage,lage,"m_age");
  M.allocate(fage,lage,"M");
  sex.allocate(fage,lage,"sex");
  prop_bay.allocate(1,tstep,fage+1,lage,"prop_bay");
  prop_coast.allocate(1,tstep,fage+1,lage,"prop_coast");
  log_sd_bay.allocate(1,tstep,fage+1,lage,"log_sd_bay");
  log_sd_coast.allocate(1,tstep,fage+1,lage,"log_sd_coast");
  acoustic_prop.allocate(1,tstep,fage+1,lage,"acoustic_prop");
  acoustic_prop_sd.allocate(1,tstep,fage+1,lage,"acoustic_prop_sd");
  use_aco_prop.allocate("use_aco_prop");
  use_age_err.allocate("use_age_err");
  age_err_a.allocate(fage,lage,fage,lage,"age_err_a");
  age_err_b.allocate(sfage_b,slage_b,sfage_b,slage_b,"age_err_b");
  age_err_c.allocate(sfage_c,slage_c,sfage_c,slage_c,"age_err_c");
  age_err_a_id.allocate(fage,lage,fage,lage,"age_err_a_id");
  age_err_b_id.allocate(sfage_b,slage_b,sfage_b,slage_b,"age_err_b_id");
  age_err_c_id.allocate(sfage_c,slage_c,sfage_c,slage_c,"age_err_c_id");
  ESS_C_bay.allocate(1,2,"ESS_C_bay");
  ESS_C_ac.allocate("ESS_C_ac");
  ESS_I_cm.allocate("ESS_I_cm");
  ESS_I_de30.allocate("ESS_I_de30");
  ESS_I_ct.allocate("ESS_I_ct");
  ESS_I_md.allocate("ESS_I_md");
  ESS_I_nj.allocate("ESS_I_nj");
  ESS_I_des.allocate("ESS_I_des");
  ESS_I_ny.allocate("ESS_I_ny");
  test.allocate("test");
  C_var.allocate(1,region,1,tstep,fmyear,lmyear);
  I_var_cm.allocate(1,tstep,fmyear,lmyear);
  I_var_de30.allocate(fmyear,lmyear);
  I_var_ct.allocate(1,tstep,fmyear,lmyear);
  I_var_md.allocate(fmyear,lmyear);
  I_var_nj.allocate(fmyear,lmyear);
  I_var_des.allocate(fmyear,lmyear);
  I_var_ny.allocate(fmyear,lmyear);
  I_var_age1.allocate(1,age1surv,fmyear,lmyear);
  I_var_age1n.allocate(1,age1surv,fmyear,lmyear);
  I_var_age1m.allocate(fmyear,lmyear);
  I_var_yoy_coast.allocate(1,yoysurv_coast,fmyear,lmyear);
  I_var_yoy_bay.allocate(1,yoysurv_bay,fmyear,lmyear);
  mod_age_err_a.allocate(sfage_a,slage_a,sfage_a,slage_a);
  mod_age_err_b.allocate(sfage_b,slage_b,sfage_b,slage_b);
  mod_age_err_c.allocate(sfage_c,slage_c,sfage_c,slage_c);
  if(test!=12345)
  {
    //Write error message and output data to the screen
    cout << "Error in model inputs!" << endl;
    cout << "data years" << endl << fdyear << endl << ldyear << endl;
    cout << "model years" << endl << fmyear << endl << lmyear << endl;
    cout << "ages" << endl << fage << endl << lage << endl;
    cout << "num years" << endl << nyrs << endl;
    cout << "fleet" << endl << fleet << endl;
    cout << "region" << endl << region << endl;
    cout << "tstep" << endl << tstep << endl;
    cout << "stock" << endl << stock << endl;
    cout << "nsurv" << endl << nsurv << endl;
    cout << "age1surv" << endl <<  age1surv << endl;
    //cout << "nage1" << endl << nage1 << endl;
    cout << "yoy surv coast"  << endl << yoysurv_coast << endl;
    cout << "yoysurv bay" << endl << yoysurv_bay << endl;
    cout << "tblock" << endl << tblock << endl;
    cout << "first tblock yr" << endl << ftyear << endl;
    //cout << "lbound of log fs" << endl << lbound << endl;
    //cout << "ubound of log fs" << endl <<  ubound  << endl;
    cout << "fage a"  << endl << sfage_a << endl;
    cout << "fage b" << endl << sfage_b << endl;
    cout << "fage c" << endl << sfage_c << endl;
    cout << "lage a"  << endl <<  slage_a << endl;
    cout << "lage b" << endl << slage_b << endl;
    cout << "lage c" << endl << slage_c << endl;
    cout << "obs catch" << endl << obs_C << endl;
    cout << "obs prop at age in catch" << endl << obs_Cp << endl;
    cout << "obs cv" << endl << obs_CV << endl;
    cout << "obs I md" << endl << obs_I_md << endl;
    cout << "obs I ny" << endl << obs_I_ny << endl;
    cout << "obs I nj" << endl << obs_I_nj << endl;
    cout << "obs I dessn" << endl << obs_I_des << endl;
    cout << "obs I de30" << endl << obs_I_de30 << endl;
    cout << "obs I cm" << endl << obs_I_cm << endl;
    cout << "obs I ct" << endl << obs_I_ct << endl;
    cout << "obs prop md" << endl << obs_Ip_md << endl;
    cout << "obs prop ny" << endl << obs_Ip_ny << endl;
    cout << "obs prop nj" << endl << obs_Ip_nj << endl;
    cout << "obs prop dessn" << endl << obs_Ip_des << endl;
    cout << "obs prop de30" << endl <<obs_Ip_de30 << endl;
    cout << "opbs prop cm" << endl <<obs_Ip_cm << endl;
    cout << "obs prop ct" << endl << obs_Ip_ct << endl;
    cout << "obs cv cm" << endl << obs_I_CV_cm << endl;
    cout << "obs cv de30" << endl << obs_I_CV_de30 << endl;
    cout << "obs cv ct" << endl <<  obs_I_CV_ct << endl;
    cout << "obs cv mdssn"  << endl << obs_I_CV_md << endl;
    cout << "obs cv njbt" << endl << obs_I_CV_nj << endl;
    cout << "obs cv dessn" << endl << obs_I_CV_des << endl;
    cout << "obs i cv ny" << endl << obs_I_CV_ny << endl;
    cout << " obs age 1 ny" << endl <<obs_I_age1n << endl;
    cout << " obs age 1 md " << endl <<  obs_I_age1m << endl;
    cout << "obs age 1ny cv " << obs_I_age1n_CV << endl;
    cout << "obs I age 1 my cv" << endl << obs_I_age1m_CV << endl;
    cout << "obs I yoy coast " << endl <<  obs_I_yoy_coast << endl;
    cout << "obs I yoy bay " << endl <<  obs_I_yoy_bay << endl;
    cout << "obs I yoy cv coast" << endl << obs_I_yoy_CV_coast << endl;
    cout << "obs I yoy cv " << endl <<  obs_I_yoy_CV_bay << endl;
    cout << "avg weight" <<  w_age << endl;
    cout << "ssbw" << endl <<ssbw_age << endl; //adjustment of rivard weight to match the time of spawning
    cout << "rivard weight" << endl << rw_age << endl; //rivard weight at age
    cout << "mat at age" << endl <<  m_age << endl;
    cout << "M " << endl << M << endl;
    cout << "sex prop" << endl <<  sex << endl;
    cout << "prop cb" << endl << prop_bay << endl;
    cout << "prop ac" << endl << prop_coast << endl;
    cout << "acoustic prop" << endl << acoustic_prop << endl;
    cout << "acoustic prop sd" << endl << acoustic_prop_sd << endl;
    cout << "age err a" << endl << age_err_a << endl;
    cout << "age err b" << endl << age_err_b << endl;
    cout << "age error c" << endl <<  age_err_c << endl;
    //cout << "ESS_C" << endl << ESS_C << endl;
    cout << "test" << endl << test << endl;
    exit(1);  //exit the program
  }
  //cout << I_var << endl;
  //convert SDs to variances so that it only has to be done once
  for(r=1;r<=region;r++)
  {
    for(t=1;t<=tstep;t++)
    {
      //cout << r << endl;
      //cout << t << endl;
      //cout << fmyear << endl;
      //cout << lmyear << endl;
      //cout << C_var(r,t) << endl;
      C_var(r,t)=square(obs_CV(r,t)(fmyear,lmyear));//calculate variance from standard deviations
      //cout << obs_CV(r,t)(fmyear,lmyear) << endl;
      //cout << C_var(r,t) << endl;
    }//close tstep
  }//close region 
  //cout << "end C_var loop"  << endl;
  
  //convert spatial SD to variance
  I_var_md=square(obs_I_CV_md(fmyear,lmyear));
  I_var_de30=square(obs_I_CV_de30(fmyear,lmyear));
  I_var_des=square(obs_I_CV_des(fmyear,lmyear));
  I_var_nj=square(obs_I_CV_nj(fmyear,lmyear));
  I_var_ny=square(obs_I_CV_ny(fmyear,lmyear));
  for(t=1;t<=tstep;t++)
  {
    I_var_ct(t)=square(obs_I_CV_ct(t)(fmyear,lmyear));
    I_var_cm(t)=square(obs_I_CV_cm(t)(fmyear,lmyear));
  }
  //cout << "end Ivar loop" << endl;
  
  for(o=1;o<=age1surv;o++)
  {
    I_var_age1n(o)=square(obs_I_age1n_CV(o)(fmyear,lmyear));
    //I_var_age1n(o)=square(obs_I_age1_CV(o)(fmyear,lmyear));
  }
  //I_var_age1n=square(obs_I_age1n_CV(fmyear,lmyear));
  //cout << "I var age1n" << endl; 
  I_var_age1m=square(obs_I_age1m_CV(fmyear,lmyear));
  //cout << "I var age1m" << endl; 
  for(z=1;z<=yoysurv_coast;z++)
  {
    I_var_yoy_coast(z)=square(obs_I_yoy_CV_coast(z)(fmyear,lmyear));
  }
    for(z=1;z<=yoysurv_bay;z++)
  {
   I_var_yoy_bay(z)=square(obs_I_yoy_CV_bay(z)(fmyear,lmyear));
  }
  //cout << "i var yoy" << endl;
  if(use_age_err==1)
  {
    mod_age_err_a=age_err_a;
    mod_age_err_b=age_err_b;
    mod_age_err_c=age_err_c;
  }
  else
  {    
    mod_age_err_a=age_err_a_id;
    mod_age_err_b=age_err_b_id;
    mod_age_err_c=age_err_c_id;
  }
  //cout << "aging error a" << endl << mod_age_err_a << endl;
  //cout << "aging error b" << endl << mod_age_err_b << endl;
  //cout << "aging error c" << endl << mod_age_err_c << endl;
  //exit(1);
  //cout << "finishvar calc " << endl;
  /*
  //cout << I_var << endl;
  cout << "data years" << endl << fdyear << endl << ldyear << endl;
  cout << "model years" << endl << fmyear << endl << lmyear << endl;
  cout << "ages" << endl << fage << endl << lage << endl;
  cout << "obs catch" << endl << obs_C << endl;
  cout << "obs prop at age in catch" << endl << obs_Cp << endl;
  cout << "obs md index" << endl << obs_I_md << endl;
  cout << "obs cm  index" << endl << obs_I_cm << endl;
  cout << "obs md prop at age in survey" << endl << obs_Ip_md << endl;
  cout << "obs cm prop at age in survey" << endl << obs_Ip_cm << endl;
  exit(1);
  */
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  log_N0_devs.allocate(fage+1,lage,-10,10,-1,"log_N0_devs");
  log_sf1_ac.allocate(1,tblock,1,tstep,-2,2,2,"log_sf1_ac");
  sf1_ac.allocate(1,tblock,1,tstep,"sf1_ac");
  #ifndef NO_AD_INITIALIZE
    sf1_ac.initialize();
  #endif
  log_sf2_ac.allocate(1,tblock,1,tstep,0,5,2,"log_sf2_ac");
  sf2_ac.allocate(1,tblock,1,tstep,"sf2_ac");
  #ifndef NO_AD_INITIALIZE
    sf2_ac.initialize();
  #endif
  log_sf1_cb.allocate(1,tblock,1,tstep,-2,2,2,"log_sf1_cb");
  sf1_cb.allocate(1,tblock,1,tstep,"sf1_cb");
  #ifndef NO_AD_INITIALIZE
    sf1_cb.initialize();
  #endif
  log_sf2_cb.allocate(1,tblock,1,tstep,0,5,2,"log_sf2_cb");
  sf2_cb.allocate(1,tblock,1,tstep,"sf2_cb");
  #ifndef NO_AD_INITIALIZE
    sf2_cb.initialize();
  #endif
  log_sf3_cb.allocate(1,tblock,1,tstep,-2,2,2,"log_sf3_cb");
  sf3_cb.allocate(1,tblock,1,tstep,"sf3_cb");
  #ifndef NO_AD_INITIALIZE
    sf3_cb.initialize();
  #endif
  log_sf4_cb.allocate(1,tblock,1,tstep,0.5,5,2,"log_sf4_cb");
  sf4_cb.allocate(1,tblock,1,tstep,"sf4_cb");
  #ifndef NO_AD_INITIALIZE
    sf4_cb.initialize();
  #endif
  log_ssf1_md.allocate(-2,2,3,"log_ssf1_md");
  ssf1_md.allocate("ssf1_md");
  #ifndef NO_AD_INITIALIZE
  ssf1_md.initialize();
  #endif
  log_ssf2_md.allocate(0,5,3,"log_ssf2_md");
  ssf2_md.allocate("ssf2_md");
  #ifndef NO_AD_INITIALIZE
  ssf2_md.initialize();
  #endif
  log_ssf1_des.allocate(-2,2,3,"log_ssf1_des");
  ssf1_des.allocate("ssf1_des");
  #ifndef NO_AD_INITIALIZE
  ssf1_des.initialize();
  #endif
  log_ssf2_des.allocate(0,5,3,"log_ssf2_des");
  ssf2_des.allocate("ssf2_des");
  #ifndef NO_AD_INITIALIZE
  ssf2_des.initialize();
  #endif
  log_ssf1_cm.allocate(1,tstep,-2,2,3,"log_ssf1_cm");
  ssf1_cm.allocate(1,tstep,"ssf1_cm");
  #ifndef NO_AD_INITIALIZE
    ssf1_cm.initialize();
  #endif
  log_ssf2_cm.allocate(1,tstep,0,5,3,"log_ssf2_cm");
  ssf2_cm.allocate(1,tstep,"ssf2_cm");
  #ifndef NO_AD_INITIALIZE
    ssf2_cm.initialize();
  #endif
  log_ssf1_ny.allocate(-2,2,3,"log_ssf1_ny");
  ssf1_ny.allocate("ssf1_ny");
  #ifndef NO_AD_INITIALIZE
  ssf1_ny.initialize();
  #endif
  log_ssf2_ny.allocate(0,5,3,"log_ssf2_ny");
  ssf2_ny.allocate("ssf2_ny");
  #ifndef NO_AD_INITIALIZE
  ssf2_ny.initialize();
  #endif
  log_ssf3_ny.allocate(-2,2,3,"log_ssf3_ny");
  ssf3_ny.allocate("ssf3_ny");
  #ifndef NO_AD_INITIALIZE
  ssf3_ny.initialize();
  #endif
  log_ssf4_ny.allocate(1,5,3,"log_ssf4_ny");
  ssf4_ny.allocate("ssf4_ny");
  #ifndef NO_AD_INITIALIZE
  ssf4_ny.initialize();
  #endif
  log_ssf1_nj.allocate(-2,2,3,"log_ssf1_nj");
  ssf1_nj.allocate("ssf1_nj");
  #ifndef NO_AD_INITIALIZE
  ssf1_nj.initialize();
  #endif
  log_ssf2_nj.allocate(0,5,3,"log_ssf2_nj");
  ssf2_nj.allocate("ssf2_nj");
  #ifndef NO_AD_INITIALIZE
  ssf2_nj.initialize();
  #endif
  log_ssf3_nj.allocate(-2,2,3,"log_ssf3_nj");
  ssf3_nj.allocate("ssf3_nj");
  #ifndef NO_AD_INITIALIZE
  ssf3_nj.initialize();
  #endif
  log_ssf4_nj.allocate(1,5,3,"log_ssf4_nj");
  ssf4_nj.allocate("ssf4_nj");
  #ifndef NO_AD_INITIALIZE
  ssf4_nj.initialize();
  #endif
  log_ssf1_de30.allocate(-2,2,3,"log_ssf1_de30");
  ssf1_de30.allocate("ssf1_de30");
  #ifndef NO_AD_INITIALIZE
  ssf1_de30.initialize();
  #endif
  log_ssf2_de30.allocate(0,5,3,"log_ssf2_de30");
  ssf2_de30.allocate("ssf2_de30");
  #ifndef NO_AD_INITIALIZE
  ssf2_de30.initialize();
  #endif
  log_ssf3_de30.allocate(-4,2,3,"log_ssf3_de30");
  ssf3_de30.allocate("ssf3_de30");
  #ifndef NO_AD_INITIALIZE
  ssf3_de30.initialize();
  #endif
  log_ssf4_de30.allocate(1,5,3,"log_ssf4_de30");
  ssf4_de30.allocate("ssf4_de30");
  #ifndef NO_AD_INITIALIZE
  ssf4_de30.initialize();
  #endif
  log_ssf1_ct.allocate(1,tstep,-2,2,3,"log_ssf1_ct");
  ssf1_ct.allocate(1,tstep,"ssf1_ct");
  #ifndef NO_AD_INITIALIZE
    ssf1_ct.initialize();
  #endif
  log_ssf2_ct.allocate(1,tstep,0,5,3,"log_ssf2_ct");
  ssf2_ct.allocate(1,tstep,"ssf2_ct");
  #ifndef NO_AD_INITIALIZE
    ssf2_ct.initialize();
  #endif
  log_ssf3_ct.allocate(1,tstep,-2,2,3,"log_ssf3_ct");
  ssf3_ct.allocate(1,tstep,"ssf3_ct");
  #ifndef NO_AD_INITIALIZE
    ssf3_ct.initialize();
  #endif
  log_ssf4_ct.allocate(1,tstep,1,5,3,"log_ssf4_ct");
  ssf4_ct.allocate(1,tstep,"ssf4_ct");
  #ifndef NO_AD_INITIALIZE
    ssf4_ct.initialize();
  #endif
  log_a_sf1.allocate(-2,2,2,"log_a_sf1");
  log_a_sf2.allocate(0,5,2,"log_a_sf2");
  a_sf1.allocate("a_sf1");
  #ifndef NO_AD_INITIALIZE
  a_sf1.initialize();
  #endif
  a_sf2.allocate("a_sf2");
  #ifndef NO_AD_INITIALIZE
  a_sf2.initialize();
  #endif
  log_q_md.allocate(1,"log_q_md");
  log_q_cm.allocate(1,tstep,1,"log_q_cm");
  log_q_ct.allocate(1,tstep,1,"log_q_ct");
  log_q_ny.allocate(1,"log_q_ny");
  log_q_nj.allocate(1,"log_q_nj");
  log_q_des.allocate(1,"log_q_des");
  log_q_de30.allocate(1,"log_q_de30");
  q_md.allocate("q_md");
  #ifndef NO_AD_INITIALIZE
  q_md.initialize();
  #endif
  q_cm.allocate("q_cm");
  q_ct.allocate("q_ct");
  q_ny.allocate("q_ny");
  #ifndef NO_AD_INITIALIZE
  q_ny.initialize();
  #endif
  q_nj.allocate("q_nj");
  #ifndef NO_AD_INITIALIZE
  q_nj.initialize();
  #endif
  q_des.allocate("q_des");
  #ifndef NO_AD_INITIALIZE
  q_des.initialize();
  #endif
  q_de30.allocate("q_de30");
  #ifndef NO_AD_INITIALIZE
  q_de30.initialize();
  #endif
  log_q_age1n.allocate(1,"log_q_age1n");
  log_q_age1m.allocate(1,"log_q_age1m");
  log_q_yoy_coast.allocate(1,yoysurv_coast,1,"log_q_yoy_coast");
  log_q_yoy_bay.allocate(1,yoysurv_bay,1,"log_q_yoy_bay");
  log_R.allocate(1,region,0,20,1,"log_R");
  log_Rdevs1.allocate(fmyear,lmyear,-10,10,1,"log_Rdevs1");
  log_Rdevs2.allocate(fmyear,lmyear,-10,10,1,"log_Rdevs2");
  log_Feq.allocate(1,stock,-5,2,2,"log_Feq");
  log_F.allocate(1,region,1,tstep,-5,2,2,"log_F");
  log_Fdevs_r1t1.allocate(fmyear,lmyear,-15,15,2,"log_Fdevs_r1t1");
  log_Fdevs_r1t2.allocate(fmyear,lmyear,-15,15,2,"log_Fdevs_r1t2");
  log_Fdevs_r2t1.allocate(fmyear,lmyear,-15,15,2,"log_Fdevs_r2t1");
  log_Fdevs_r2t2.allocate(fmyear,lmyear,-15,15,2,"log_Fdevs_r2t2");
  Feq.allocate(1,stock,fage,lage,"Feq");
  #ifndef NO_AD_INITIALIZE
    Feq.initialize();
  #endif
  Nt.allocate(1,stock,1,tstep,fmyear,lmyear,fage,lage,"Nt");
  #ifndef NO_AD_INITIALIZE
    Nt.initialize();
  #endif
  N.allocate(1,stock,1,region,1,tstep,fmyear,lmyear,fage,lage,"N");
  #ifndef NO_AD_INITIALIZE
    N.initialize();
  #endif
  N_region.allocate(1,region,1,tstep,fmyear,lmyear,fage,lage,"N_region");
  #ifndef NO_AD_INITIALIZE
    N_region.initialize();
  #endif
  Nbay.allocate(1,tstep,fmyear,lmyear,fage,lage,"Nbay");
  #ifndef NO_AD_INITIALIZE
    Nbay.initialize();
  #endif
  Ncoast.allocate(1,tstep,fmyear,lmyear,fage,lage,"Ncoast");
  #ifndef NO_AD_INITIALIZE
    Ncoast.initialize();
  #endif
  Z.allocate(1,region,1,tstep,fmyear,lmyear,fage,lage,"Z");
  #ifndef NO_AD_INITIALIZE
    Z.initialize();
  #endif
  F.allocate(1,region,1,tstep,fmyear,lmyear,fage,lage,"F");
  #ifndef NO_AD_INITIALIZE
    F.initialize();
  #endif
  Fbar.allocate(1,stock,fmyear,lmyear,fage,lage,"Fbar");
  #ifndef NO_AD_INITIALIZE
    Fbar.initialize();
  #endif
  Fplus.allocate(1,stock,fmyear,lmyear,"Fplus");
  #ifndef NO_AD_INITIALIZE
    Fplus.initialize();
  #endif
  Freg_7plus.allocate(1,region,1,tstep,fmyear,lmyear,7,lage,"Freg_7plus");
  #ifndef NO_AD_INITIALIZE
    Freg_7plus.initialize();
  #endif
  fsel.allocate(1,region,1,tstep,fmyear,lmyear,fage,lage,"fsel");
  #ifndef NO_AD_INITIALIZE
    fsel.initialize();
  #endif
  log_fsel.allocate(1,region,1,tstep,fmyear,lmyear,fage,lage,"log_fsel");
  #ifndef NO_AD_INITIALIZE
    log_fsel.initialize();
  #endif
  afsel.allocate(fage,lage,"afsel");
  #ifndef NO_AD_INITIALIZE
    afsel.initialize();
  #endif
  log_afsel.allocate(fage,lage,"log_afsel");
  #ifndef NO_AD_INITIALIZE
    log_afsel.initialize();
  #endif
  ssel_md.allocate(sfage_b,slage_b,"ssel_md");
  #ifndef NO_AD_INITIALIZE
    ssel_md.initialize();
  #endif
  ssel_cm.allocate(1,tstep,sfage_a,slage_a,"ssel_cm");
  #ifndef NO_AD_INITIALIZE
    ssel_cm.initialize();
  #endif
  ssel_ct.allocate(1,tstep,sfage_a,slage_a,"ssel_ct");
  #ifndef NO_AD_INITIALIZE
    ssel_ct.initialize();
  #endif
  ssel_ny.allocate(sfage_c,slage_c,"ssel_ny");
  #ifndef NO_AD_INITIALIZE
    ssel_ny.initialize();
  #endif
  ssel_nj.allocate(sfage_b,slage_b,"ssel_nj");
  #ifndef NO_AD_INITIALIZE
    ssel_nj.initialize();
  #endif
  ssel_des.allocate(sfage_c,slage_c,"ssel_des");
  #ifndef NO_AD_INITIALIZE
    ssel_des.initialize();
  #endif
  ssel_de30.allocate(sfage_a,slage_a,"ssel_de30");
  #ifndef NO_AD_INITIALIZE
    ssel_de30.initialize();
  #endif
  q_age1n.allocate("q_age1n");
  #ifndef NO_AD_INITIALIZE
  q_age1n.initialize();
  #endif
  q_age1m.allocate("q_age1m");
  #ifndef NO_AD_INITIALIZE
  q_age1m.initialize();
  #endif
  q_yoy_coast.allocate(1,yoysurv_coast,"q_yoy_coast");
  #ifndef NO_AD_INITIALIZE
    q_yoy_coast.initialize();
  #endif
  q_yoy_bay.allocate(1,yoysurv_bay,"q_yoy_bay");
  #ifndef NO_AD_INITIALIZE
    q_yoy_bay.initialize();
  #endif
  q_age1.allocate(1,age1surv,"q_age1");
  #ifndef NO_AD_INITIALIZE
    q_age1.initialize();
  #endif
  est_C.allocate(1,stock,1,region,1,tstep,fmyear,lmyear,"est_C");
  #ifndef NO_AD_INITIALIZE
    est_C.initialize();
  #endif
  est_region_C.allocate(1,region,1,tstep,fmyear,lmyear,"est_region_C");
  #ifndef NO_AD_INITIALIZE
    est_region_C.initialize();
  #endif
  est_C_age.allocate(1,stock,1,region,1,tstep,fmyear,lmyear,fage,lage,"est_C_age");
  #ifndef NO_AD_INITIALIZE
    est_C_age.initialize();
  #endif
  est_C_age_err.allocate(1,stock,1,region,1,tstep,fmyear,lmyear,fage,lage,"est_C_age_err");
  #ifndef NO_AD_INITIALIZE
    est_C_age_err.initialize();
  #endif
  est_totC_age.allocate(1,region,1,tstep,fmyear,lmyear,fage,lage,"est_totC_age");
  #ifndef NO_AD_INITIALIZE
    est_totC_age.initialize();
  #endif
  est_totC_age_err.allocate(1,region,1,tstep,fmyear,lmyear,fage,lage,"est_totC_age_err");
  #ifndef NO_AD_INITIALIZE
    est_totC_age_err.initialize();
  #endif
  est_Cp.allocate(1,stock,1,region,1,tstep,fmyear,lmyear,fage,lage,"est_Cp");
  #ifndef NO_AD_INITIALIZE
    est_Cp.initialize();
  #endif
  est_region_Cp.allocate(1,region,1,tstep,fmyear,lmyear,fage,lage,"est_region_Cp");
  #ifndef NO_AD_INITIALIZE
    est_region_Cp.initialize();
  #endif
  sigma2_Cp.allocate(1,region,1,tstep,fmyear,lmyear,fage,lage,"sigma2_Cp");
  #ifndef NO_AD_INITIALIZE
    sigma2_Cp.initialize();
  #endif
  ESS_C.allocate(1,region,fmyear,lmyear,"ESS_C");
  #ifndef NO_AD_INITIALIZE
    ESS_C.initialize();
  #endif
  est_I_age_md.allocate(fmyear,lmyear,sfage_b,slage_b,"est_I_age_md");
  #ifndef NO_AD_INITIALIZE
    est_I_age_md.initialize();
  #endif
  est_I_age_md_err.allocate(fmyear,lmyear,sfage_b,slage_b,"est_I_age_md_err");
  #ifndef NO_AD_INITIALIZE
    est_I_age_md_err.initialize();
  #endif
  est_Ip_md.allocate(fmyear,lmyear,sfage_b,slage_b,"est_Ip_md");
  #ifndef NO_AD_INITIALIZE
    est_Ip_md.initialize();
  #endif
  est_I_md.allocate(fmyear,lmyear,"est_I_md");
  #ifndef NO_AD_INITIALIZE
    est_I_md.initialize();
  #endif
  sigma2_Ip_md.allocate(fmyear,lmyear,sfage_b,slage_b,"sigma2_Ip_md");
  #ifndef NO_AD_INITIALIZE
    sigma2_Ip_md.initialize();
  #endif
  est_I_age_cm.allocate(1,tstep,fmyear,lmyear,sfage_a,slage_a,"est_I_age_cm");
  #ifndef NO_AD_INITIALIZE
    est_I_age_cm.initialize();
  #endif
  est_Ip_cm.allocate(1,tstep,fmyear,lmyear,sfage_a,slage_a,"est_Ip_cm");
  #ifndef NO_AD_INITIALIZE
    est_Ip_cm.initialize();
  #endif
  est_I_cm.allocate(1,tstep,fmyear,lmyear,"est_I_cm");
  #ifndef NO_AD_INITIALIZE
    est_I_cm.initialize();
  #endif
  sigma2_Ip_cm.allocate(1,tstep,fmyear,lmyear,sfage_a,slage_a,"sigma2_Ip_cm");
  #ifndef NO_AD_INITIALIZE
    sigma2_Ip_cm.initialize();
  #endif
  est_I_age_ct.allocate(1,tstep,fmyear,lmyear,sfage_a,slage_a,"est_I_age_ct");
  #ifndef NO_AD_INITIALIZE
    est_I_age_ct.initialize();
  #endif
  est_I_age_ct_err.allocate(1,tstep,fmyear,lmyear,sfage_a,slage_a,"est_I_age_ct_err");
  #ifndef NO_AD_INITIALIZE
    est_I_age_ct_err.initialize();
  #endif
  est_Ip_ct.allocate(1,tstep,fmyear,lmyear,sfage_a,slage_a,"est_Ip_ct");
  #ifndef NO_AD_INITIALIZE
    est_Ip_ct.initialize();
  #endif
  est_I_ct.allocate(1,tstep,fmyear,lmyear,"est_I_ct");
  #ifndef NO_AD_INITIALIZE
    est_I_ct.initialize();
  #endif
  sigma2_Ip_ct.allocate(1,tstep,fmyear,lmyear,sfage_a,slage_a,"sigma2_Ip_ct");
  #ifndef NO_AD_INITIALIZE
    sigma2_Ip_ct.initialize();
  #endif
  est_I_age_ny.allocate(fmyear,lmyear,sfage_c,slage_c,"est_I_age_ny");
  #ifndef NO_AD_INITIALIZE
    est_I_age_ny.initialize();
  #endif
  est_I_age_ny_err.allocate(fmyear,lmyear,sfage_c,slage_c,"est_I_age_ny_err");
  #ifndef NO_AD_INITIALIZE
    est_I_age_ny_err.initialize();
  #endif
  est_Ip_ny.allocate(fmyear,lmyear,sfage_c,slage_c,"est_Ip_ny");
  #ifndef NO_AD_INITIALIZE
    est_Ip_ny.initialize();
  #endif
  est_I_ny.allocate(fmyear,lmyear,"est_I_ny");
  #ifndef NO_AD_INITIALIZE
    est_I_ny.initialize();
  #endif
  sigma2_Ip_ny.allocate(fmyear,lmyear,sfage_c,slage_c,"sigma2_Ip_ny");
  #ifndef NO_AD_INITIALIZE
    sigma2_Ip_ny.initialize();
  #endif
  est_I_age_nj.allocate(fmyear,lmyear,sfage_b,slage_b,"est_I_age_nj");
  #ifndef NO_AD_INITIALIZE
    est_I_age_nj.initialize();
  #endif
  est_I_age_nj_err.allocate(fmyear,lmyear,sfage_b,slage_b,"est_I_age_nj_err");
  #ifndef NO_AD_INITIALIZE
    est_I_age_nj_err.initialize();
  #endif
  est_Ip_nj.allocate(fmyear,lmyear,sfage_b,slage_b,"est_Ip_nj");
  #ifndef NO_AD_INITIALIZE
    est_Ip_nj.initialize();
  #endif
  est_I_nj.allocate(fmyear,lmyear,"est_I_nj");
  #ifndef NO_AD_INITIALIZE
    est_I_nj.initialize();
  #endif
  sigma2_Ip_nj.allocate(fmyear,lmyear,sfage_b,slage_b,"sigma2_Ip_nj");
  #ifndef NO_AD_INITIALIZE
    sigma2_Ip_nj.initialize();
  #endif
  est_I_age_des.allocate(fmyear,lmyear,sfage_c,slage_c,"est_I_age_des");
  #ifndef NO_AD_INITIALIZE
    est_I_age_des.initialize();
  #endif
  est_I_age_des_err.allocate(fmyear,lmyear,sfage_c,slage_c,"est_I_age_des_err");
  #ifndef NO_AD_INITIALIZE
    est_I_age_des_err.initialize();
  #endif
  est_Ip_des.allocate(fmyear,lmyear,sfage_c,slage_c,"est_Ip_des");
  #ifndef NO_AD_INITIALIZE
    est_Ip_des.initialize();
  #endif
  est_I_des.allocate(fmyear,lmyear,"est_I_des");
  #ifndef NO_AD_INITIALIZE
    est_I_des.initialize();
  #endif
  sigma2_Ip_des.allocate(fmyear,lmyear,sfage_c,slage_c,"sigma2_Ip_des");
  #ifndef NO_AD_INITIALIZE
    sigma2_Ip_des.initialize();
  #endif
  est_I_age_de30.allocate(fmyear,lmyear,sfage_a,slage_a,"est_I_age_de30");
  #ifndef NO_AD_INITIALIZE
    est_I_age_de30.initialize();
  #endif
  est_Ip_de30.allocate(fmyear,lmyear,sfage_a,slage_a,"est_Ip_de30");
  #ifndef NO_AD_INITIALIZE
    est_Ip_de30.initialize();
  #endif
  est_I_de30.allocate(fmyear,lmyear,"est_I_de30");
  #ifndef NO_AD_INITIALIZE
    est_I_de30.initialize();
  #endif
  sigma2_Ip_de30.allocate(fmyear,lmyear,sfage_a,slage_a,"sigma2_Ip_de30");
  #ifndef NO_AD_INITIALIZE
    sigma2_Ip_de30.initialize();
  #endif
  est_I_age1n.allocate(1,age1surv,fmyear,lmyear,"est_I_age1n");
  #ifndef NO_AD_INITIALIZE
    est_I_age1n.initialize();
  #endif
  est_I_age1m.allocate(fmyear,lmyear,"est_I_age1m");
  #ifndef NO_AD_INITIALIZE
    est_I_age1m.initialize();
  #endif
  est_I_yoy_coast.allocate(1,yoysurv_coast,fmyear,lmyear,"est_I_yoy_coast");
  #ifndef NO_AD_INITIALIZE
    est_I_yoy_coast.initialize();
  #endif
  est_I_yoy_bay.allocate(1,yoysurv_bay,fmyear,lmyear,"est_I_yoy_bay");
  #ifndef NO_AD_INITIALIZE
    est_I_yoy_bay.initialize();
  #endif
  est_I_age1.allocate(1,age1surv,fmyear,lmyear,"est_I_age1");
  #ifndef NO_AD_INITIALIZE
    est_I_age1.initialize();
  #endif
  rw.allocate(fmyear,lmyear,fage,lage,"rw");
  #ifndef NO_AD_INITIALIZE
    rw.initialize();
  #endif
  SSB_w.allocate(fmyear,lmyear,fage,lage,"SSB_w");
  #ifndef NO_AD_INITIALIZE
    SSB_w.initialize();
  #endif
  B.allocate(1,stock,1,region,fmyear,lmyear,"B");
  #ifndef NO_AD_INITIALIZE
    B.initialize();
  #endif
  SSB.allocate(1,stock,1,region,fmyear,lmyear,"SSB");
  #ifndef NO_AD_INITIALIZE
    SSB.initialize();
  #endif
  J_w.allocate(fmyear,lmyear,"J_w");
  #ifndef NO_AD_INITIALIZE
    J_w.initialize();
  #endif
  log_prop_bay.allocate(1,tstep,fage+1,lage,-10,0,5,"log_prop_bay");
  log_prop_coast.allocate(1,tstep,fage+1,lage,-10,0,5,"log_prop_coast");
  prop.allocate(1,stock,1,tstep,1,region,fage,lage,"prop");
  #ifndef NO_AD_INITIALIZE
    prop.initialize();
  #endif
  Ntot_cb.allocate(fmyear,lmyear,"Ntot_cb");
  Ntot_ac.allocate(fmyear,lmyear,"Ntot_ac");
  N_cbincb1.allocate(fmyear,lmyear,"N_cbincb1");
  N_cbincb2.allocate(fmyear,lmyear,"N_cbincb2");
  N_cbinac1.allocate(fmyear,lmyear,"N_cbinac1");
  N_cbinac2.allocate(fmyear,lmyear,"N_cbinac2");
  N_acincb1.allocate(fmyear,lmyear,"N_acincb1");
  N_acincb2.allocate(fmyear,lmyear,"N_acincb2");
  N_acinac1.allocate(fmyear,lmyear,"N_acinac1");
  N_acinac2.allocate(fmyear,lmyear,"N_acinac2");
  logSSB_cb.allocate(fmyear,lmyear,"logSSB_cb");
  logSSB_ac.allocate(fmyear,lmyear,"logSSB_ac");
  recruit_cb.allocate(fmyear,lmyear,"recruit_cb");
  recruit_ac.allocate(fmyear,lmyear,"recruit_ac");
  Fstockcb_sd.allocate(fmyear,lmyear,"Fstockcb_sd");
  Fstockac_sd.allocate(fmyear,lmyear,"Fstockac_sd");
  Fcb1_sd.allocate(fmyear,lmyear,"Fcb1_sd");
  Fcb2_sd.allocate(fmyear,lmyear,"Fcb2_sd");
  Fac1_sd.allocate(fmyear,lmyear,"Fac1_sd");
  Fac2_sd.allocate(fmyear,lmyear,"Fac2_sd");
  mdssn.allocate(fmyear,lmyear,"mdssn");
  chesmmap1.allocate(fmyear,lmyear,"chesmmap1");
  chesmmap2.allocate(fmyear,lmyear,"chesmmap2");
  njbt.allocate(fmyear,lmyear,"njbt");
  nyohs.allocate(fmyear,lmyear,"nyohs");
  dessn.allocate(fmyear,lmyear,"dessn");
  de30.allocate(fmyear,lmyear,"de30");
  ctlist1.allocate(fmyear,lmyear,"ctlist1");
  ctlist2.allocate(fmyear,lmyear,"ctlist2");
  oceanage1.allocate(fmyear,lmyear,"oceanage1");
  bayage1.allocate(fmyear,lmyear,"bayage1");
  cbyoy.allocate(fmyear,lmyear-1,"cbyoy");
  nyyoy.allocate(fmyear,lmyear-1,"nyyoy");
  njyoy.allocate(fmyear,lmyear-1,"njyoy");
  sig2_f.allocate(1,region,1,tstep,"sig2_f");
  #ifndef NO_AD_INITIALIZE
    sig2_f.initialize();
  #endif
  Lcatch.allocate(1,region,1,tstep,"Lcatch");
  #ifndef NO_AD_INITIALIZE
    Lcatch.initialize();
  #endif
  Lcatchagecomp.allocate(1,region,1,tstep,"Lcatchagecomp");
  #ifndef NO_AD_INITIALIZE
    Lcatchagecomp.initialize();
  #endif
  Lindex_md.allocate("Lindex_md");
  #ifndef NO_AD_INITIALIZE
  Lindex_md.initialize();
  #endif
  Lindexagecomp_md.allocate("Lindexagecomp_md");
  #ifndef NO_AD_INITIALIZE
  Lindexagecomp_md.initialize();
  #endif
  Lindex_ny.allocate("Lindex_ny");
  #ifndef NO_AD_INITIALIZE
  Lindex_ny.initialize();
  #endif
  Lindexagecomp_ny.allocate("Lindexagecomp_ny");
  #ifndef NO_AD_INITIALIZE
  Lindexagecomp_ny.initialize();
  #endif
  Lindex_nj.allocate("Lindex_nj");
  #ifndef NO_AD_INITIALIZE
  Lindex_nj.initialize();
  #endif
  Lindexagecomp_nj.allocate("Lindexagecomp_nj");
  #ifndef NO_AD_INITIALIZE
  Lindexagecomp_nj.initialize();
  #endif
  Lindex_des.allocate("Lindex_des");
  #ifndef NO_AD_INITIALIZE
  Lindex_des.initialize();
  #endif
  Lindexagecomp_des.allocate("Lindexagecomp_des");
  #ifndef NO_AD_INITIALIZE
  Lindexagecomp_des.initialize();
  #endif
  Lindex_de30.allocate("Lindex_de30");
  #ifndef NO_AD_INITIALIZE
  Lindex_de30.initialize();
  #endif
  Lindexagecomp_de30.allocate("Lindexagecomp_de30");
  #ifndef NO_AD_INITIALIZE
  Lindexagecomp_de30.initialize();
  #endif
  Lindex_ct.allocate(1,tstep,"Lindex_ct");
  #ifndef NO_AD_INITIALIZE
    Lindex_ct.initialize();
  #endif
  Lindexagecomp_ct.allocate(1,tstep,"Lindexagecomp_ct");
  #ifndef NO_AD_INITIALIZE
    Lindexagecomp_ct.initialize();
  #endif
  Lindex_cm.allocate(1,tstep,"Lindex_cm");
  #ifndef NO_AD_INITIALIZE
    Lindex_cm.initialize();
  #endif
  Lindexagecomp_cm.allocate(1,tstep,"Lindexagecomp_cm");
  #ifndef NO_AD_INITIALIZE
    Lindexagecomp_cm.initialize();
  #endif
  Lage1index.allocate(1,age1surv,"Lage1index");
  #ifndef NO_AD_INITIALIZE
    Lage1index.initialize();
  #endif
  Lage1index_ny.allocate(1,age1surv,"Lage1index_ny");
  #ifndef NO_AD_INITIALIZE
    Lage1index_ny.initialize();
  #endif
  Lage1index_md.allocate("Lage1index_md");
  #ifndef NO_AD_INITIALIZE
  Lage1index_md.initialize();
  #endif
  Lyoyindex_coast.allocate(1,yoysurv_coast,"Lyoyindex_coast");
  #ifndef NO_AD_INITIALIZE
    Lyoyindex_coast.initialize();
  #endif
  Lyoyindex_bay.allocate(1,yoysurv_bay,"Lyoyindex_bay");
  #ifndef NO_AD_INITIALIZE
    Lyoyindex_bay.initialize();
  #endif
  Lfsel.allocate(1,region,1,tstep,"Lfsel");
  #ifndef NO_AD_INITIALIZE
    Lfsel.initialize();
  #endif
  Lssel_md.allocate("Lssel_md");
  #ifndef NO_AD_INITIALIZE
  Lssel_md.initialize();
  #endif
  Lssel_cm.allocate(1,tstep,"Lssel_cm");
  #ifndef NO_AD_INITIALIZE
    Lssel_cm.initialize();
  #endif
  Lssel_ny.allocate("Lssel_ny");
  #ifndef NO_AD_INITIALIZE
  Lssel_ny.initialize();
  #endif
  Lssel_nj.allocate("Lssel_nj");
  #ifndef NO_AD_INITIALIZE
  Lssel_nj.initialize();
  #endif
  Lssel_ct.allocate(1,tstep,"Lssel_ct");
  #ifndef NO_AD_INITIALIZE
    Lssel_ct.initialize();
  #endif
  Lssel_des.allocate("Lssel_des");
  #ifndef NO_AD_INITIALIZE
  Lssel_des.initialize();
  #endif
  Lssel_de30.allocate("Lssel_de30");
  #ifndef NO_AD_INITIALIZE
  Lssel_de30.initialize();
  #endif
  Lssel_a.allocate("Lssel_a");
  #ifndef NO_AD_INITIALIZE
  Lssel_a.initialize();
  #endif
  Lssel_b.allocate("Lssel_b");
  #ifndef NO_AD_INITIALIZE
  Lssel_b.initialize();
  #endif
  Lssel_c.allocate("Lssel_c");
  #ifndef NO_AD_INITIALIZE
  Lssel_c.initialize();
  #endif
  pen_N0_dev.allocate("pen_N0_dev");
  #ifndef NO_AD_INITIALIZE
  pen_N0_dev.initialize();
  #endif
  pen_f2sel.allocate("pen_f2sel");
  #ifndef NO_AD_INITIALIZE
  pen_f2sel.initialize();
  #endif
  pen_F.allocate("pen_F");
  #ifndef NO_AD_INITIALIZE
  pen_F.initialize();
  #endif
  pen_prop.allocate("pen_prop");
  #ifndef NO_AD_INITIALIZE
  pen_prop.initialize();
  #endif
  pen_prop_aco.allocate("pen_prop_aco");
  #ifndef NO_AD_INITIALIZE
  pen_prop_aco.initialize();
  #endif
  pen_prop_bay.allocate("pen_prop_bay");
  #ifndef NO_AD_INITIALIZE
  pen_prop_bay.initialize();
  #endif
  pen_cb_sel.allocate("pen_cb_sel");
  #ifndef NO_AD_INITIALIZE
  pen_cb_sel.initialize();
  #endif
  pen_rdev.allocate("pen_rdev");
  #ifndef NO_AD_INITIALIZE
  pen_rdev.initialize();
  #endif
  pen_fdev.allocate("pen_fdev");
  #ifndef NO_AD_INITIALIZE
  pen_fdev.initialize();
  #endif
  pen_feq.allocate("pen_feq");
  #ifndef NO_AD_INITIALIZE
  pen_feq.initialize();
  #endif
  pen_sf.allocate("pen_sf");
  #ifndef NO_AD_INITIALIZE
  pen_sf.initialize();
  #endif
  pen_sf_ct.allocate("pen_sf_ct");
  #ifndef NO_AD_INITIALIZE
  pen_sf_ct.initialize();
  #endif
  neg_LL.allocate("neg_LL");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  //set starting values for parameters
  //log_N0=log(200000.);//mean abudance
  log_R=log(1000000.); //mean of log recruitment
  log_F(1)=log(0.3); //mean of log F
  log_F(2)=log(0.3); //mean of log F
  log_Feq=log(0.3); //F in the first year
  //*******************
  //starting parameters for selectivity
  //*******************
  //fishing selectivity
  log_afsel=0.;  //average fsel
  //atlantic coast
  log_sf1_ac(1)=0.;
  log_sf1_ac(2)=0.;
  log_sf2_ac(1)=2.;
  log_sf2_ac(2)=2.;
  //chesapeake bay
  log_sf1_cb(1)=0.;
  log_sf1_cb(2)=0.;
  log_sf2_cb(1)=2.;
  log_sf2_cb(2)=2.;
  log_sf3_cb(1)=0.;
  log_sf3_cb(2)=0.;
  log_sf4_cb(1)=2.7;
  log_sf4_cb(2)=2.7;
  //log_fs_cb(1)=0.; //fsel for ches bay time step 1
  //log_fs_cb(2)=0.; //fsel for ches bay time step 2
  //log_fs_ac(1)=0.; //fsel for atl coast time step 1
  //log_fs_ac(2)=0.; //fsel for atl coast time step 2
  //survey selectivity
  //mdssn
  log_ssf1_md=0.;
  log_ssf2_md=2.;
  //chemmap
  log_ssf1_cm(1)=0.;
  log_ssf1_cm(2)=0.;
  log_ssf2_cm(1)=1;
  log_ssf2_cm(2)=1;
  //log_ssf3_cm(1)=0.;
  //log_ssf3_cm(2)=0.;
  //log_ssf4_cm(1)=2.3;
  //log_ssf4_cm(2)=2.3;
  //dessn
  log_ssf1_des=0.;
  log_ssf2_des=2.;
  //ctlist
  log_ssf1_ct(1)=0.;
  log_ssf1_ct(2)=0.;
  log_ssf2_ct(1)=1.;
  log_ssf2_ct(2)=1.;
  log_ssf3_ct(1)=0.;
  log_ssf3_ct(2)=0.;
  log_ssf4_ct(1)=2.7;
  log_ssf4_ct(2)=2.7;
  //nyohs
  log_ssf1_ny=0.;
  log_ssf2_ny=2.;
  log_ssf3_ny=0.;
  log_ssf4_ny=2.5;
  //njbt
  log_ssf1_nj=0.;
  log_ssf2_nj=2.;
  log_ssf3_nj=0.;
  log_ssf4_nj=2.5;
  //de30
  log_ssf1_de30=0.;
  log_ssf2_de30=2.;
  log_ssf3_de30=0.;
  log_ssf4_de30=2.5;
  
  log_a_sf1=log(1.);
  log_a_sf2=log(4.);
  log_q_md=log(0.000005); //log of md SSN survey catchability
  log_q_cm=log(0.000005); //log of chesmmap surv catchability
  log_q_ct=log(0.000005); //log of CT LIST surv catchability
  log_q_ny=log(0.000005); // log of NY OHS catchability
  log_q_nj=log(0.000005); //log of nybt catchability
  log_q_des=log(0.000005); // log of del SSN, catchability
  log_q_de30=log(0.000005); //log of de30 catchability
  log_q_age1n=log(0.000005);       //using this until age 1 ny is disaggregated
  log_q_age1m=log(0.000005);       //log MD age 1 survey catchability
  log_q_yoy_coast=log(0.000005);   //log YOY catchability for coast (NY and NJ agg)
  log_q_yoy_bay=log(0.000005);    //log YOY catchability for bay (MD and VA)
  //defining log_prop_bay and coast so that you are not taking the log of a 0, instead value equal -10
  for(ts=1;ts<=tstep;ts++)
  {
    for(a=fage+1;a<=lage;a++)
    {
      log_prop_bay(ts,a)=prop_bay(ts,a);
    }//close age loop
  }//close timestep loop
  for(ts=1;ts<=tstep;ts++)
  {
    for(a=fage+1;a<=lage;a++)
    {
      log_prop_coast(ts,a)=prop_coast(ts,a);
    }//close age loop
  }//close tstep loop
  for(y=fmyear;y<=lmyear;y++)
  {
    r=1;
    if(y<1990)
      {
        ESS_C(r,y)=ESS_C_bay(1); //setting ESS_C to a differnt ESS during Moratorium in the bay
      }//close if statement
      else
      {
        ESS_C(r,y)=ESS_C_bay(2);
      }//close else statement
     r=2;
     ESS_C(r,y)=ESS_C_ac;
   }
  //cout << "prop bay in coast" << endl << exp(log_prop_bay) << endl;
  //cout << "prop coast in coast" << endl << exp(log_prop_coast) << endl;
  //cout << age_err_a << endl;
  //exit(1);
  //cout << "finish ss and q calc " << endl;
  
}

void model_parameters::userfunction(void)
{
  neg_LL =0.0;
  calculate_mortality();
  //cout << "After calculate_mortality" <<endl;
  calculate_N_at_age();
  //cout << "After calculate_N_at_age" <<endl;
  calculate_F_stock();
  //cout << "After calculate_F_stock" <<endl;
  calculate_fishery_catch();
  //cout << "After calculate_fisher_catch" <<endl;
  calculate_survey_catch();
  //cout << "After calculate_survey_catch" <<endl;
  evaluate_likelihood();  
  //cout << "After calculate_likelihood" <<endl;
  calculate_B_SSB();
  //cout << "After B_SSB" << endl;
}

void model_parameters::calculate_mortality(void)
{
  //*******************************
  // Calculate Fishery Selectivity
  //******************************
  sf1_ac=mfexp(log_sf1_ac);
  sf2_ac=mfexp(log_sf2_ac);
  sf1_cb=mfexp(log_sf1_cb);
  sf2_cb=mfexp(log_sf2_cb);
  sf3_cb=mfexp(log_sf3_cb);
  sf4_cb=mfexp(log_sf4_cb);
  for(t=1;t<=tblock;t++)
  {
    ftbyr=ftyear(t);
    if(t==tblock)
    {
      ftlyr=lmyear;
    }
    else
    {
      ftlyr=ftyear(t+1)-1;
    }
    //cout << "First year timeblock" << ftbyr << endl;
    //cout << "Last Year timeblock" << ftlyr << endl;
  /*
  // arc tangent calcualtion to bound fisheries selectivity between -10, 1
   for(y=ftbyr;y<=ftlyr;y++)
    {
      for(ts=1;ts<=tstep;ts++)
      {
      for(a=fage;a<=lage-1;a++)
      {
        atan_log_fs_cb(t,ts,a)=atan(log_fs_cb(t,ts,a));//arc tangent transformation of fsel, this transforms values to -pi/2 to pi/2
        atan_log_fs_cb(t,ts,a)=atan_log_fs_cb(t,ts,a)/(PI/2);//dividing arctan transformation to make between -1 and 1
        atan_log_fs_cb(t,ts,a)=atan_log_fs_cb(t,ts,a)*((ubound-lbound)/2)+(lbound/2);//this widens the range to be 12 and shifts negative to be between -10 and 2
        atan_log_fs_ac(t,ts,a)=atan(log_fs_ac(t,ts,a));//arc tangent transformation of fsel, this transforms values to -pi/2 to pi/2
        atan_log_fs_ac(t,ts,a)=atan_log_fs_ac(t,ts,a)/(PI/2);//dividing arctan transformation to make between -1 and 1
        atan_log_fs_ac(t,ts,a)=atan_log_fs_ac(t,ts,a)*((ubound-lbound)/2)+(lbound/2);//this widens the range to be 12 and shifts negative to be between -10 and 2
     }//close age loop
    }//close timeblock loop
  }//close fleet loop
  */
  //cout << "end atan calc" << endl;
  //cout << "fsel" << " " << fsel << endl; 
    r=1;//fill fsel for chesapeake bay
    for(y=ftbyr;y<=ftlyr;y++)
    {
      for(ts=1;ts<=tstep;ts++)
      {
        /*
        for(a=fage;a<=6;a++)
        {
          //cout << y << " " << " yr" <<  endl;
          //cout << r << " " << " reg" << endl;
          //cout << a << " " << " age" << endl;
          //cout << t << " " << "tbloc" << endl;
          //cout << ts << " " << "tstep" << endl;
          //cout << atan_log_fs(r,t,a) << endl;
          fsel(r,ts,y,a)=exp(atan_log_fs_cb(t,ts,a));
          //cout << fsel(y,ts,a,r) << endl;
        }//close age loop for fsel
        fsel(r,ts,y,a)=1.0; //estimating selectivity relative to age 7
        for(a=8;a<=lage;a++)
        {
          //cout << a << endl; 
          fsel(r,ts,y,a)=exp(atan_log_fs_cb(t,ts,a-1));
         // cout << fsel(y,ts,a,r) << endl;
        }//close sel after age 8 loop
        */
        for(a=fage;a<=lage;a++)
        {
          fsel(r,ts,y,a)=(1./(1.+mfexp(-sf1_cb(t,ts)*(double(a)-sf2_cb(t,ts)))))*
                  (1./(1.+mfexp(-sf3_cb(t,ts)*(sf4_cb(t,ts)- double(a)))));
        }//close age loop
        fsel(r,ts,y)/=max(fsel(r,ts,y));
      }//close tstep
    } //close year loop for fsel
    r=2;//fill fsel for atlantic coast
    for(y=ftbyr;y<=ftlyr;y++)
    {
      for(ts=1;ts<=tstep;ts++)
      {
        /*
        for(a=fage;a<=6;a++)
        {
          //cout << y << " " << " yr" <<  endl;
          //cout << r << " " << " reg" << endl;
          //cout << a << " " << " age" << endl;
          //cout << t << " " << "tbloc" << endl;
          //cout << ts << " " << "tstep" << endl;
          //cout << atan_log_fs(r,t,a) << endl;
          fsel(r,ts,y,a)=exp(atan_log_fs_ac(t,ts,a));
          //cout << fsel(y,ts,a,r) << endl;
        }//close age loop for fsel
        fsel(r,ts,y,a)=1.0; //estimating selectivity relative to age 7
        for(a=8;a<=lage;a++)
        {
          //cout << a << endl; 
          fsel(r,ts,y,a)=exp(atan_log_fs_ac(t,ts,a-1));
         // cout << fsel(y,ts,a,r) << endl;
        }//close sel after age 8 loop
        */
        for(a=fage;a<=lage;a++)
        {
          fsel(r,ts,y,a)=1./(1.+mfexp(-sf1_ac(t,ts)*(double(a)-sf2_ac(t,ts))));
        }//close age loop
        fsel(r,ts,y)/=max(fsel(r,ts,y));
      }//close tstep
   } //close year loop for fsel
  }//close time block
  //cout << fsel << endl;
  //
  //exit(1);
  //cout << "fsel" << endl;
  for(y=fmyear;y<=lmyear;y++)
  {
    for(r=1;r<=region;r++)
    {
      for(ts=1;ts<=tstep;ts++)
      {
        log_fsel(r,ts,y)=log(fsel(r,ts,y));
      }//clost ts loop
    }//close region loop
  }//close year loop
  //cout << log_fsel << endl;
  //exit(1);
  //cout << "fsel" << fsel << endl;
  //exit(1);
  //cout << "end fsel" << endl;
  //*******************************
  // Calculate Survey Selectivity
  //******************************
  //SN added logistic or double logistic 5/25/23 to restrict model
  //group a - ages 1 - 15
  //set starting parameters
  ssf1_des=mfexp(log_ssf1_des);
  ssf2_des=mfexp(log_ssf2_des);
  ssf1_md=mfexp(log_ssf1_md);
  ssf2_md=mfexp(log_ssf2_md);
  ssf1_cm=mfexp(log_ssf1_cm);
  ssf2_cm=mfexp(log_ssf2_cm);
  //ssf3_cm=mfexp(log_ssf3_cm);
  //ssf4_cm=mfexp(log_ssf4_cm);
  ssf1_ct=mfexp(log_ssf1_ct);
  ssf2_ct=mfexp(log_ssf2_ct);
  ssf3_ct=mfexp(log_ssf3_ct);
  ssf4_ct=mfexp(log_ssf4_ct);
  ssf1_ny=mfexp(log_ssf1_ny);
  ssf2_ny=mfexp(log_ssf2_ny);
  ssf3_ny=mfexp(log_ssf3_ny);
  ssf4_ny=mfexp(log_ssf4_ny);
  ssf1_nj=mfexp(log_ssf1_nj);
  ssf2_nj=mfexp(log_ssf2_nj);
  ssf3_nj=mfexp(log_ssf3_nj);
  ssf4_nj=mfexp(log_ssf4_nj);
  ssf1_de30=mfexp(log_ssf1_de30);
  ssf2_de30=mfexp(log_ssf2_de30);
  ssf3_de30=mfexp(log_ssf3_de30);
  ssf4_de30=mfexp(log_ssf4_de30);
  //cout << "end starting params" << endl;
  //calculating ssel
  for(a=sfage_a;a<=slage_a;a++)
  {
    //de30 , double logisitic
    ssel_de30(a)=(1./(1.+mfexp(-ssf1_de30*(double(a)-ssf2_de30))))*
                  (1./(1.+mfexp(-ssf3_de30*(ssf4_de30-double(a)))));
  //cout << ssel_de30 << endl << "end ssel de30" << endl;
  //cout << ssel_cm << endl;
      for(t=1;t<=tstep;t++)
      {
      //ct list, double logistic
      ssel_ct(t,a)=(1./(1.+mfexp(-ssf1_ct(t)*(double(a)-ssf2_ct(t)))))*
      (1./(1.+mfexp(-ssf3_ct(t)*(ssf4_ct(t)-double(a)))));
      //chesmmap,  logistic
      ssel_cm(t,a)=(1./(1.+mfexp(-ssf1_cm(t)*(double(a)-ssf2_cm(t)))));//*
                  //(1./(1.+mfexp(-ssf3_cm(t)*(ssf4_cm(t)-double(a)))));//turned off double logisitic 11/9/12
      //cout << "age" << " " << a << endl << "ts" << " " << ts << endl;
      }//close tstep
  }//close age loop
  //cout << ssel_ct << endl;
  ssel_de30/=max(ssel_de30);
  for(t=1;t<=tstep;t++)
  {
    ssel_ct(t)/=max(ssel_ct(t));
    ssel_cm(t)/=max(ssel_cm(t));
  }
  //cout << "end ssel group a" << endl;
  //group b ages 2-15
  for(a=sfage_b;a<=slage_b;a++)
  {
    //md ssn, logistic
    ssel_md(a)=1./(1.+mfexp(-ssf1_md*(double(a)-ssf2_md)));
    //njbt, double logistic
    ssel_nj(a)=(1./(1.+mfexp(-ssf1_nj*(double(a)-ssf2_nj))))*
                  (1./(1.+mfexp(-ssf3_nj*(ssf4_nj- double(a)))));
  }//close age loop
  ssel_md/=max(ssel_md);
  ssel_nj/=max(ssel_nj);
  //cout << "end ssel group b" << endl;
  //group c, age 2-13
  for(a=sfage_c;a<=slage_c;a++)
  {
    //nyohs, double logistic
    ssel_ny(a)=(1./(1.+mfexp(-ssf1_ny*(double(a)-ssf2_ny))))*
                  (1./(1.+mfexp(-ssf3_ny*(ssf4_ny- double(a)))));
    //dessn, logisitic
    ssel_des(a)=1./(1.+mfexp(-ssf1_des*(double(a)-ssf2_des)));
  }//close age loop
  ssel_ny/=max(ssel_ny);
  ssel_des/=max(ssel_des);
  //cout << "end ssel group c" << endl;
  //exit(1);
  /*
  //no function, selectivity is free for all ages
  for(a=sfage_a;a<=6;a++)
  {
    //cout << "a" << " " << a << endl;
    ssel_de30(a)=exp(log_ss_de30(a));
    //cout << "de30" << " " << ssel_de30 << endl;
    for(t=1;t<=tstep;t++)
    {
     //cout << "t" << " " << t << endl;
      ssel_cm(t,a)=exp(log_ss_cm(t,a));
      ssel_ct(t,a)=exp(log_ss_ct(t,a));
     }//close age loop for timesetp
   }//close age loop
   ssel_de30(7)=1.0; //estimating selectivity relative to age 7
   for(t=1;t<=tstep;t++)
   {
     ssel_ct(t,7)=1.0;
     ssel_cm(t,7)=1.0;
   }//close time step loop
     for(a=8;a<=slage_a;a++)
     {
      ssel_de30(a)=exp(log_ss_de30(a-1));
      for(t=1;t<=tstep;t++)
      {
        ssel_ct(t,a)=exp(log_ss_ct(t,a-1));
        ssel_cm(t,a)=exp(log_ss_cm(t,a-1));
      }//close tstep loop
   }//close age loop
  //cout << "end cm and dt ssel" << endl;
  for(a=sfage_b;a<=6;a++)
  {
    ssel_md(a)=exp(log_ss_md(a));
    ssel_nj(a)=exp(log_ss_nj(a));
  }//close age loop for fsel
  ssel_md(7)=1.0; //estimating selectivity relative to age 7
  ssel_nj(7)=1.0; 
  for(a=8;a<=slage_b;a++)
  {
    ssel_md(a)=exp(log_ss_md(a-1));
    ssel_nj(a)=exp(log_ss_nj(a-1));
  }//close age loop
  //  }//close survey loop
  //cout << ssel_md << endl;
  //exit(1);
  //Group C ssel
  //for(s=1;s<=nsurv_c;s++)//looping over surveys and age
  //{
  for(a=sfage_c;a<=6;a++)
  {
    ssel_ny(a)=exp(log_ss_ny(a));
    ssel_des(a)=exp(log_ss_des(a));
  }//close age loop
  ssel_ny(7)=1.0; //estimating selectivity relative to age 7
  ssel_des(7)=1.0;
  for(a=8;a<=slage_c;a++)
  {
    ssel_ny(a)=exp(log_ss_ny(a-1));
    ssel_des(a)=exp(log_ss_des(a-1));
  }//close age loop
  // }//close survey loop
  */
  //cout << "end ssel" << endl;
  //Catchability
  q_md=mfexp(log_q_md); //log of md SSN survey catchability
  q_cm=mfexp(log_q_cm); //log of chesmmap surv catchability
  q_ct=mfexp(log_q_ct); //log of CT LIST surv catchability
  q_ny=mfexp(log_q_ny); // log of NY OHS catchability
  q_nj=mfexp(log_q_nj); //log of nybt catchability
  q_des=mfexp(log_q_des); // log of del SSN, catchability
  q_de30=mfexp(log_q_de30); //log of de30 catchability
  q_age1n=mfexp(log_q_age1n);
  q_age1m=mfexp(log_q_age1m);
  for(z=1;z<=yoysurv_coast;z++)
  {
    q_yoy_coast(z)=mfexp(log_q_yoy_coast(z));
  }
  for(z=1;z<=yoysurv_bay;z++)
  {
    q_yoy_bay(z)=mfexp(log_q_yoy_bay(z));  
  }
  /*
  cout << q_a << endl;
  cout << q_b << endl;
  cout << q_c << endl;
  exit(1);
  */
  //calculate the fishing mortality rate (F) for each year and age
  a_sf1=mfexp(log_a_sf1);
  a_sf2=mfexp(log_a_sf2);
  for(a=fage;a<=lage;a++)
  {
    //a_sf1 = slope; a_sf2=age at 50% selc
    //afsel(a)=1./(1.+mfexp(-1.*(double(a)-4.)));//setting afsel to a logistic curve with slope of 1 and 50% sel at 4
    afsel(a)=1./(1.+mfexp(-a_sf1*(double(a)-a_sf2)));//setting afsel to a logistic curve with slope of 1 and 50% sel at 4
  }
  //afsel(fage,lage)/=max(afsel);
  //afsel=(fsel(1,1,fmyear)+fsel(1,2,fmyear)+fsel(2,1,fmyear)+fsel(2,2,fmyear))/4.0;//average selectivity in the first year
  //cout << "a-sf1= " << a_sf1 << endl << "a_sf2= " << a_sf2 << endl;
  //cout << "afsel= " << afsel << endl;
  //exit(1);
  //fishing eq in the first year
  for(s=1;s<=stock;s++)
  {
    Feq(s)=afsel*mfexp(log_Feq(s));
  }//close stock loop
  //cout << Feq << endl;
  //exit(1);
  //calculate F in the first year, second time step for region 1
  r=1;
  t=2;
  for(a=fage;a<=lage;a++)
  {
     //F(r,t,fmyear,a)=fsel(r,t,fmyear,a)*exp(log_F(r,t)+log_Fdevs_r1t2(fmyear));
  }//close age
  //Fill out F for the rest of years for time step 1 region 1
  t=1;
  for(y=fmyear;y<=lmyear;y++)
  {       
    for(a=fage;a<=lage;a++)
    {
        F(r,t,y,a)=fsel(r,t,y,a)*mfexp(log_F(r,t)+log_Fdevs_r1t1(y));
    }//close age loop
  }//close year loop
  //cout << F <<  endl;
  //exit(1);
  //fill out rest of years for  
  t=2;
  for(y=fmyear;y<=lmyear;y++)
  {       
    for(a=fage;a<=lage;a++)
    {
      F(r,t,y,a)=fsel(r,t,y,a)*mfexp(log_F(r,t)+log_Fdevs_r1t2(y));
    }//close age loop
  }//close year loop
  //fill out first year  of F for region 2, timestep 2
  r=2;
  t=2;
  for(a=fage;a<=lage;a++)
  {
     //F(r,t,fmyear,a)=fsel(r,t,fmyear,a)*exp(log_F(r,t)+log_Fdevs_r2t2(fmyear));
    }//close age
  //Fill out F for the rest of years fore region 2, time step 1
  t=1;
  for(y=fmyear;y<=lmyear;y++)
  {       
    for(a=fage;a<=lage;a++)
    {
      F(r,t,y,a)=fsel(r,t,y,a)*mfexp(log_F(r,t)+log_Fdevs_r2t1(y));
    }//close age loop
  }//close year loop
  //Fill out F for the rest of years fore region 2, time step 2
  t=2;
  for(y=fmyear;y<=lmyear;y++)
  {       
    for(a=fage;a<=lage;a++)
    {
      F(r,t,y,a)=fsel(r,t,y,a)*mfexp(log_F(r,t)+log_Fdevs_r2t2(y));
    }//close age loop
  }//close year loop
  //cout << F << endl;
  //exit(1);
  /*
  for(t=1;t<=tstep;t++)
    {
      for(r=1;r<=region;r++)
      {
          cout << "t=" << t << ",r=" << r << endl << F(r,t) << endl;
      }
   }
  */
  //calculate the total mortality rate for each year and age 
  for(y=fmyear;y<=lmyear;y++)
  {
    for(t=1;t<=tstep;t++)
    {
      for(r=1;r<=region;r++)
      {
        Z(r,t,y)=M; //M divided by 2 in dat file to account for time step
        Z(r,t,y)+=F(r,t,y);
      }//close region
    }//close tstep
  }//close year
  //cout << "Z" << endl << Z << endl << endl;
  //cout << "end mort" << endl;
  //exit(1);
}

void model_parameters::calculate_N_at_age(void)
{
  //filling in the prop matrix
  /*
  for(ts=1;ts<=tstep;ts++)
  {
    for(a=fage;a<=lage;a++)
    {
      prop(1,ts,1,a)=exp(log_prop_bay(ts,a));//chesapeake bay stock in cb
      prop(2,ts,1,a)=exp(log_prop_coast(ts,a));//atlantic coast stock in cb
    }
  }
  for(ts=1;ts<=tstep;ts++)
  {
    for(a=fage;a<=lage;a++)
    {
      prop(1,ts,2,a)=1.-exp(log_prop_bay(ts,a));//chesapeake bay stock in cb
      prop(2,ts,2,a)=1.-exp(log_prop_coast(ts,a));//atlantic coast stock in cb
    }
  }
 */
  //filling occupancy probabilities for ages 1 and 2 so that fish cannot migrate at these ages
  for(ts=1;ts<=tstep;ts++)
  { 
    prop(1,ts,1,1)=1.;//100% of bay stock at age 1 in bay
    prop(1,ts,2,1)=0.;//0% of bay stock in coast at age 1
    prop(2,ts,1,1)=0.;//0% of coast stock in bay at age 1
    prop(2,ts,2,1)=1.;//100% of coast stock in coast at age 1
  }
  //cout << "prop" << endl << prop << endl;
  //exit(1);
  for(s=1;s<=stock;s++)
  {
    for(ts=1;ts<=tstep;ts++)
    {
      for(r=1;r<=region;r++)
      {
        for(a=fage+1;a<=lage;a++)
        {
          if(r==1)//if region = chesapeake bay
          {
            prop(1,ts,r,a)=1.-mfexp(log_prop_bay(ts,a));//chesapeake bay stock in cb
            prop(2,ts,r,a)=1.-mfexp(log_prop_coast(ts,a));//atlantic coast stock in cb
            //prop(s,ts,r,a)=exp(atan(log_prop_bay(ts,a))/(PI/2)*((0-(-10))/2)+(-10/2));// arc tangent transformation of occupancy probabilities to have bounded between log(-10) and log(0)
          }
          else
          {
            prop(1,ts,r,a)=mfexp(log_prop_bay(ts,a));//chesapeake bay stock in cb
            prop(2,ts,r,a)=mfexp(log_prop_coast(ts,a));//atlantic coast stock in cb
          } //end else statement
        }//close age loop
      }//close region loop
    }//close ts loop
  }//close stock loop
  //cout << "prop" << endl << prop << endl;
  //exit(1);
  //fill in abundance at age in the first row of the N-at-age matrix
  //fill in recuitment for the first year
  Nt(1,1,fmyear,fage)=mfexp(log_R(1)+log_Rdevs1(fmyear));//Bay stock (stock, tstep, year,,age); mfexp allows exp of moderate numbers to be differentiable
  Nt(2,1,fmyear,fage)=mfexp(log_R(2)+log_Rdevs2(fmyear));//Atlantic stock(stock, tstep, year,,age)  //*MW* using the same Rdevs each year.
  //cout << "recruit" << endl;
  //calculate N at age in the first year and first time step, assuming that population is at equillibrium
  for(s=1;s<=stock;s++)
  {
    for(a=fage+1;a<=lage;a++)
    {
      //cout << "a" << " " << a << endl;
      Nt(s,1,fmyear,a)=Nt(s,1,fmyear,a-1)*mfexp(-(2.*M(a-1)+Feq(s,a-1)));//M*2 because we are getting ages from age 1 to age2, not from one time step to the next +F(2,fmyear,a-1)*fsel(2,fmyear,a-1)));
      //cout << Nt(s,1,fmyear,a) << endl;
    }//close age loop
    //including plus group in first year
    Nt(s,1,fmyear,lage)/=1.-mfexp(-(2.*M(lage)+Feq(s,lage)));//+F(2,fmyear,lage)*fsel(2,fmyear,lage)));
    //now include N0 deviations in first year of equilibrium
    Nt(s,1,fmyear)(fage+1,lage)=elem_prod(Nt(s,1,fmyear)(fage+1,lage),mfexp(log_N0_devs(fage+1,lage)));
    //*MW* modified
    for(r=1;r<=region;r++)
    {
      //calculate abundance in each region in the first year time step 1
      N(s,r,1,fmyear)=migrate(Nt(s,1,fmyear),prop(s,1,r)); //*MW* 
    }//close region loop
  }//close stock loop
  //cout << afsel << endl;
  //exit(1);
  //cout << Nt << endl;
  /*
  for(s=1;s<=stock;s++)
  {
    for(r=1;r<=region;r++)
    {
      for(t=1;t<=2;t++)
      {
        cout << "s=" << s << ", r=" << r << ", t=" << t << endl;
        cout << N(s,r,t,fmyear) << endl;
      }
    }
  }
  exit(1);
  */
  //cout << "end N tstep 1 calc first year" << endl;
  //Fill in abundance in first year, time step 2 *MW* modified below.
  for(s=1;s<=stock;s++)
  {
    for(a=fage;a<=lage;a++)
    { 
      //Nt(s,2,fmyear,a)=(N(s,1,1,fmyear,a)*mfexp(-(M(a)+F(1,1,fmyear,a)*fsel(1,1,fmyear,a)))+N(s,2,1,fmyear,a)*mfexp(-(M(a)+F(2,1,fmyear,a)*fsel(2,1,fmyear,a))));
      Nt(s,2,fmyear,a)=(N(s,1,1,fmyear,a)*mfexp(-(M(a)+F(1,1,fmyear,a)))+N(s,2,1,fmyear,a)*mfexp(-(M(a)+F(2,1,fmyear,a))));//SN edited 2/19/24 because were applying fsel twice
    }//close age loop
    for(r=1;r<=region;r++)
    {
      N(s,r,2,fmyear)=migrate(Nt(s,2,fmyear),prop(s,2,r)); //
    }//close region loop
  }//close stock loop
  /*
  cout << Nt << endl << endl;
  for(s=1;s<=stock;s++)
  {
    for(r=1;r<=region;r++)
    {
      for(t=1;t<=2;t++)
      {
        cout << "s=" << s << ", r=" << r << ", t=" << t << endl;
        cout << N(s,r,t) << endl;
      }
    }
  }
  exit(1);
  */
  //cout << "end N time step 2 fmyear" << endl;
  //fill in the rest of the rows of the N-at-age matrix by applying exponential mortality
  //*MW* need to keep track of which time variable is incrementing for the different seasons
  //*MW* if t=1, then the next time step is the same year, but t=2
  //*MW* if t=2, then the next time step is t=1, y+1
  for(y=fmyear+1;y<=lmyear;y++)
  {
     //calculate recruitment for first time step of each year
     Nt(1,1,y,fage)=mfexp(log_R(1)+log_Rdevs1(y));//Bay stock, (yr,tstep,stock,age)  //*MW*
     Nt(2,1,y,fage)=mfexp(log_R(2)+log_Rdevs2(y));//coast stock, (yr,tstep,stock,age)  //*MW*
     //cout << Nt(1,1,y,fage) << endl;
     //cout << Nt(2,1,y,fage) << endl;
  }
  //cout << Nt << endl;
  //exit(1);
  /*
  for(r=1;r<=region;r++)
  {
    for(t=1;t<=2;t++)
    {
      cout << "r=" << r << ", t=" << t << endl;
      cout << Z(r,t) << endl;
    }
  }
  exit(1);
  */
  for(y=fmyear+1;y<=lmyear;y++)
  {
     for(t=1;t<=tstep;t++)
     {
       if(t==1)//first dothe first tstep for each year; need to refer to the previous year and age
       {
         t2=2; //set the time step for reference on the RHS to 2 (the one before our current time step)
         for(s=1;s<=stock;s++)
         {
           for(a=fage+1;a<=lage;a++)
           {
             Nt(s,t,y,a)=N(s,1,t2,y-1,a-1)*mfexp(-Z(1,t2,y-1,a-1))+N(s,2,t2,y-1,a-1)*mfexp(-Z(2,t2,y-1,a-1));
             //cout << "s=" << s << ", t=" << t << ", y=" << y << ", a=" << a << endl;
             //cout << Nt(s,t,y,a) << " " << N(s,1,t2,y-1,a-1) << " " <<exp(-Z(1,t2,y-1,a-1)) << " " <<N(s,2,t2,y-1,a-1) << " " <<exp(-Z(2,t2,y-1,a-1))<< endl;
           }
           //plus group
           Nt(s,1,y,lage)+=N(s,1,t2,y-1,lage)*mfexp(-Z(1,t2,y-1,lage))+ N(s,2,t2,y-1,lage)*mfexp(-Z(2,t2,y-1,lage));
           for(r=1;r<=region;r++)
           {
             //cout << a << endl;
             //cout << r << endl;
             N(s,r,t,y)=migrate(Nt(s,t,y),prop(s,t,r));
             //cout << endl << "s=" << s <<", r=" << r << ", t=" << t << ", y=" << y << endl;
             //cout << "migrate" << endl << N(s,r,t,y) << endl << Nt(s,t,y) << endl << prop(s,r,t) << endl << endl;
           }//close region loop
         }//close stock loop
       }//close if statement
       else
       {
         //t2=1;//set the time step for reference on the RHS to 1 (the one before our current time step)
         for(s=1;s<=stock;s++)
         {
           Nt(s,t,y)=elem_prod(N(s,1,t-1,y),mfexp(-Z(1,t-1,y)))+ elem_prod(N(s,2,t-1,y),mfexp(-Z(2,t-1,y)));
           //cout << "s=" << s << ", t=" << t << ", y=" << y << endl;
           //cout << Nt(s,t,y) << endl << N(s,1,t-1,y) << endl << exp(-Z(1,t-1,y)) << endl <<  N(s,2,t-1,y) << endl << exp(-Z(2,t-1,y)) << endl;
           //including plus group *MW* no extra plus group calculation for this time step because the animals only increment an age between timestep 2 and 1 of the next year.
           //Nt(s,2,y,lage)+=(N(s,1,1,y,lage)*exp(-M(lage)+F(1,1,y,lage)*fsel(1,1,y,lage))+N(s,2,1,y,lage)*exp(-M(lage)+F(2,1,y,lage)*fsel(2,1,y,lage)));
           for(r=1;r<=region;r++)
           {
             //calculate abundance in each region in the first year
             N(s,r,t,y)=migrate(Nt(s,t,y),prop(s,t,r)); // propar should be read in?
           }//close region loop
         }//close stock loop
       }//close else statement
     }//close t loop
  }//close year loop
 // for(s=1;s<=stock;s++)
 // {
 //   for(t=1;t<=2;t++)
  //  {
  //    cout << "s=" << s << ", t=" << t << endl;
  //    cout << Nt(s,t) << endl << endl;
   // }
 // }
  //for(s=1;s<=stock;s++)
  //{
  /*
  for(r=1;r<=region;r++)
    {
      for(t=1;t<=2;t++)
      {
        //cout << "s=" << s << ", r=" << r << ", t=" << t << endl;
        //cout << N(s,r,t) << endl;
        cout <<  " r=" << r << ", t=" << t << endl;
        cout << Nt(r,t) << endl;
    //  }
    }
  }
  exit(1);
  */
  //cout << "N calcs" << endl;
  //exit(1);
  //calculate the population in each region at each time step first
  for(y=fmyear;y<=lmyear;y++)
  {
    for(t=1;t<=tstep;t++)
    {
      //for(a=fage;a<=lage;a++)
      //{
        N_region(1,t,y)=N(1,1,t,y)+N(2,1,t,y); //N in bay in tstep 1 and 2
        N_region(2,t,y)=N(1,2,t,y)+N(2,2,t,y); //N in coast in tstep 1 and 2
    } //close tstep loop
  }//close year loop
  /*
  if(last_phase())
  {
  for(r=1;r<=region;r++)
  {
    for(t=1;t<=tstep;t++)
    {
       cout << "t=" << t << endl << "r=" << r << endl << N_region(r,t) << endl;
    } //close tstep loop
  }//close region loop
  //exit(1);
  }
  */
  ///Fill in sd report values
  //tot pop
  for(y=fmyear;y<=lmyear;y++)
  {
    Ntot_cb(y)=log(sum(N(1,1,1,y)+N(1,2,1,y)));//cb stock in both region, at jan1
    Ntot_ac(y)=log(sum(N(2,1,1,y)+N(2,2,1,y))); //ac stock in both regions, at jan1
  }
  //pop 7+
  for(y=fmyear;y<=lmyear;y++)
  {
    N_cbincb1(y)=log(sum(N(1,1,1,y)(7,lage)));//cb stock in cb region, timestep=1
    N_cbincb2(y)=log(sum(N(1,1,2,y)(7,lage))); //cb stock in cb region, timestep=2
    N_cbinac1(y)=log(sum(N(1,2,1,y)(7,lage))); //cb stock in ocean region, timestep=1
    N_cbinac2(y)=log(sum(N(1,2,2,y)(7,lage))); //cb stock in ocean region, timestep=2
    N_acincb1(y)=log(sum(N(2,1,1,y)(7,lage))); //ac stock in cb region, timestep =1
    N_acincb2(y)=log(sum(N(2,1,2,y)(7,lage))); //ac stock in cb region, timestep =2
    N_acinac1(y)=log(sum(N(2,2,1,y)(7,lage))); //ac stock in ac region, timestep =1
    N_acinac2(y)=log(sum(N(2,2,2,y)(7,lage))); //ac stock in ac region, timestep =2
  }//close year loop
  //age1
  for(y=fmyear;y<=lmyear;y++)
  {
    recruit_cb(y)=log(N(1,1,1,y,1));
    recruit_ac(y)=log(N(2,2,1,y,1));
  }
}

void model_parameters::calculate_B_SSB(void)
{
  for(y=fmyear;y<=lmyear;y++)
  {
    for(s=1;s<=stock;s++)
    {
      for(r=1;r<=region;r++)
      {
        B(s,r,y)=N(s,r,1,y)*rw_age(y)/1000; //January 1 Biomass; using catch  weight at age
        SSB(s,r,y)=N(s,r,1,y)*elem_prod(sex,elem_prod(ssbw_age(y),m_age))/1000; //Female Spawning Stock biomass
      }//close region loop
    }//close stock loop
  }//close year loop
  //cout << "finish ssb calcs" << endl;
    // * is the dot product, which multiplies the values and then sums it
  //fill sd report params
  for(y=fmyear;y<=lmyear;y++)
  {
    logSSB_cb(y)=log(SSB(1,1,y)+SSB(1,2,y)); //cb stock in both regions
    logSSB_ac(y)=log(SSB(2,1,y)+SSB(2,2,y)); //ac stock in both region
  }//close year loop
  /*
  for(y=fmyear;y<=lmyear;y++)
  {
    logSSB_cbincb(y)=log(SSB(1,1,y)); //cb stock in CB
    logSSB_cbinac(y)=log(SSB(1,2,y)); //cb stock in ac region
      logSSB_acincb(y)=log(SSB(2,1,y)); //ac stock in cb region
    logSSB_acinac(y)=log(SSB(2,2,y)); //ac stock in ac region
  }//close year loop
  */
}

void model_parameters::calculate_fishery_catch(void)
{
  //calculate fishery catch at age using the Baranov catch equation C=(F/Z)*(1-exp(-Z))*N
  for(y=fmyear;y<=lmyear;y++)
  {
    for(t=1;t<=tstep;t++)
    {
      for(r=1;r<=region;r++)
      {
        //cout << "r" << " " << r << endl;
        //cout << "yr" << " " << y << endl;
        est_region_C(r,t,y)=0.0;
        for(s=1;s<=stock;s++)
        {
          //cout << "s" << " " << s << endl;
          est_C_age(s,r,t,y)(fage,lage)=elem_prod(elem_prod(elem_div(F(r,t,y)(fage,lage),Z(r,t,y)(fage,lage)),1.-mfexp(-Z(r,t,y)(fage,lage))),N(s,r,t,y)(fage,lage));
          //cout << est_C_age (s,r,t,y) << endl;
	  est_C(s,r,t,y)=sum(est_C_age(s,r,t,y));  //calculate total catch by stock
          est_C_age_err(s,r,t,y)=mod_age_err_a*est_C_age(s,r,t,y);
          //cout << "r=" << r << endl;
          //cout << "t=" << t << endl;
          //cout << "est_C=" << est_C(s,r,t,y) << endl;
        }//close stock loop
        est_totC_age(r,t,y)(fage,lage)=elem_prod(elem_prod(elem_div(F(r,t,y)(fage,lage),Z(r,t,y)(fage,lage)),1.-mfexp(-Z(r,t,y)(fage,lage))),N_region(r,t,y)(fage,lage));
        est_totC_age_err(r,t,y)=mod_age_err_a*est_totC_age(r,t,y);//calculations to get catch at age for corrected aging error
        est_region_C(r,t,y)+=sum(est_totC_age_err(r,t,y));  //calculate total catch by stock
                 //use est_regon_C  for the likelihood
      } //close region loop
    }//close time step
  }//close year loop for catch at age
  //cout << "est c age" << endl << est_totC_age << endl;
  //cout << "est c age err" << endl << est_totC_age_err << endl;
  //exit(1);
  /*
  for(t=1;t<=tstep;t++)
  {
    for(r=1;r<=region;r++)
    {
      cout << "t=" << t << ",r=" << r << endl << "est_totC_age" << " " << endl <<  est_totC_age(r,t) <<endl;
      cout << "est_totC_age_err" << " " << endl << est_totC_age_err(r,t) << endl;
    }
  }
  exit(1);
  */
  /*
  for(t=1;t<=tstep;t++)
  {
    for(r=1;r<=region;r++)
    {
        cout << "r=" << r << endl << "t=" << t << endl << "reg_c" << est_region_C(r,t) << endl << "region age" << est_totC_age(r,t) << endl ;
     }
  }
  for(t=1;t<=tstep;t++)
  {
   for(r=1;r<=region;r++)
  {
      for(s=1;s<=stock;s++)
     {
        cout << "r=" << r << endl << "t=" << t << endl << est_totC_age(r,t) << endl;
      }
     }
   }
  exit(1);
  */
  //cout << "fish catch " << endl;
  //calculate proportions at age in the catch
  for(y=fmyear;y<=lmyear;y++)
  { 
    for(t=1;t<=tstep;t++)
    { 
      for(r=1;r<=region;r++)
      {
        est_region_Cp(r,t,y)=0.0;
	for(s=1;s<=stock;s++)
        {
          //cout << s << endl;
	  //est_Cp(s,r,t,y)=sum(est_totC_age_err(s,r,t,y)); //sum catch ove age
          //est_Cp(s,r,t,y)(fage,lage)/=est_totC(s,r,t,y); //divide by total catch to get proportions at age
        }
        for(a=fage;a<=lage;a++)
        {
          //cout << "t=" << t << endl << "r=" << r << endl;
          est_region_Cp(r,t,y,a)+=est_totC_age_err(r,t,y,a); //calculating proportions at age for an entire region, regardless of stock
          est_region_Cp(r,t,y,a)/=est_region_C(r,t,y);  //continue calculating prop at age for an entire region, regardless of stock
          //cout << est_region_Cp(r,t,y) << endl;
          }//close age loop
      }//close region loop
    }// close time step
  }//close year loop for estcp
  /*
  for(t=1;t<=tstep;t++)
  { 
    for(r=1;r<=region;r++)
    {
      cout << "r=," << r << "t=" << t  << endl << est_region_Cp(r,t) << endl;
    }
  }
  exit(1);
  */
  //calculate variance of expected proportion for catch at age (Fournier - multifanciel)
  for(y=fmyear;y<=lmyear;y++)
  {
    for(r=1;r<=region;r++)
    {
      if(y<ftyear(2) && r==2)
      {
        d=0.05;
      }//close if statement
      else
      {
        d=0.1;
      }//close else statement
       for(t=1;t<=tstep;t++)
       { 
         for(a=fage;a<=lage;a++)
         {
           sigma2_Cp(r,t,y,a)=((1.-est_region_Cp(r,t,y,a))*est_region_Cp(r,t,y,a)+(d/double(lage-fage+1)))/ESS_C(r,y);
         }//close age loop for sigma2cp
       }//close timestep loop
     }//close region
   }//close year
  //cout << "finish C calcs" << endl;
  /*
  cout << "obs CP" << endl;
  cout << obs_Cp << endl;
  cout << "est CP" << endl;  
  cout << est_Cp << endl;
  cout << "sigmaCP" << endl;
  cout << sigma2_Cp << endl;
  exit(1);
  */
}

void model_parameters::calculate_F_stock(void)
{
  //here we are calculating the weighted F
  //Freg_7plus.initialize(); // Initialize the result array
  //for(s=1;s<=stock;s++)
  //{
    for(r=1;r<=region;r++)
    {
      for(ts=1;ts<=tstep;ts++)
      {
        for(y=fmyear;y<=lmyear;y++)
        {
          for(a=7;a<=lage;a++)
          {
            Freg_7plus(r,ts,y,a)=F(r,ts,y,a)*(N_region(r,ts,y,a))/sum(N_region(r,ts,y)(7,lage));//(sum(N(1,r,ts,y)(4,lage))+sum(N(2,r,ts,y))); //N is abundance of each stock in each region, so you want to add abundance for both stocks
          }
        }
      }//close year
    }//clsoe region
  //}
  //cout << "N1" << endl << N << endl;
  //cout << "Freg" << endl << Freg_7plus << endl;
  //exit(1);
  for(y=fmyear;y<=lmyear;y++)
  {
    Fcb1_sd(y)=log(sum(Freg_7plus(1,1,y))); //f mort for cb region, jan-jun
    Fcb2_sd(y)=log(sum(Freg_7plus(1,2,y))); //f mort for cb region, jul-dec
    Fac1_sd(y)=log(sum(Freg_7plus(2,1,y))); //f mort for ac region, jan-jun
    Fac2_sd(y)=log(sum(Freg_7plus(2,2,y))); //f mort for ac region, jul-dec
  }//close year loop for sd report
  //here we are stimating annual F for each stock
  Fbar.initialize(); // Initialize the result array
  for(s=1;s<=stock;s++)
  {
    for(y=fmyear;y<=lmyear;y++)
    {
      for(a=fage;a<=lage;a++)
      {
        for(r=1;r<=region;r++)
        {
          for(ts=1;ts<=tstep;ts++)
          {
            Fbar(s,y,a)+=prop(s,ts,r,a)*F(r,ts,y,a);          
          }
        }
      }
    }
  }
  //cout << "Fbar" << " " << Fbar << endl;
  //calculte F for each stock of fish age 7 and older
  Fplus.initialize(); // Initialize the result array
  for(s=1;s<=stock;s++)
  {
    for(y=fmyear;y<=lmyear;y++)
    {
      for(a=7;a<=lage;a++)
      {
        Fplus(s,y)+=(Fbar(s,y,a)*Nt(s,1,y,a))/sum(Nt(s,1,y)(7,lage)); //Nt(s,1,y)is the stock abundance, in timestep 1, for each year   
      }
    }
  }
  //cout << "Fplus" << endl << Fplus << endl;
  //exit(1);
  Fstockcb_sd=log(Fplus(1)); //f mort for cb stock
  Fstockac_sd=log(Fplus(2)); //f mort for ac stock
}

void model_parameters::calculate_survey_catch(void)
{
  //calculate the population in each region at each time step first
  for(y=fmyear;y<=lmyear;y++)
  {
    for(t=1;t<=tstep;t++)
    {
      //for(a=fage;a<=lage;a++)
      //{
        //Nbay(t,y)=N(1,1,t,y)+N(2,1,t,y); //N in bay in tstep 1 and 2
        //Ncoast(t,y)=N(1,2,t,y)+N(2,2,t,y); //N in coast in tstep 1 and 2
      //} // close age loop
        r=1;
        Nbay(t,y)=N_region(r,t,y); //N in bay in tstep 1 and 2
        r=2;
        Ncoast(t,y)=N_region(r,t,y); //N in coast in tstep 1 and 2
      //cout << N(1,1,t,y) << endl << endl << N(2,1,t,y) << endl << endl;
      //cout << N(1,2,t,y) << endl << endl << N(2,2,t,y) << endl << endl;
    } //close tstep loop
  }//close year loop
  //cout << Nbay << endl << endl;
  //cout << Ncoast << endl;
  //exit(1);
  //cout << "nbay ncooast" << endl;
  ////////////////////////
  //  Chesapeake bay    //
  ////////////////////////
  //MD Spawning Stock Survey
  //calculate estimated survey catch at age for each year
  for(y=fmyear;y<=lmyear;y++)
  {
    //est_I_age_md(y)=q_md*elem_prod(Nbay(1,y)(sfage_b,slage_b),ssel_md);//timstep1, region1
    est_I_age_md(y)=q_md*elem_prod(N(1,1,1,y)(sfage_b,slage_b),ssel_md);//stock1,region1timstep1, 
    est_I_age_md_err(y)=mod_age_err_b*est_I_age_md(y);// est_i_age*aging error because MD uses otolith
  }//closing year loop
  //cout << "q_MD" << q_md << "ssel_MD"<< ssel_md << endl << Nbay <<  endl;
  //calculate total survey catch at age for each year
  //cout << "est i age md" << endl;
  est_I_md=rowsum(est_I_age_md_err);
  //cout <<"est i md" << endl << est_I_age_md << endl<<  "est i md err " <<endl << est_I_age_md_err << endl;
  //cout << est_I_md << endl;
  //exit(1);
  for(y=fmyear;y<=lmyear;y++)
  {//calculate proportions at age in the survey catch
    est_Ip_md(y)=est_I_age_md_err(y)/est_I_md(y);
    for(a=sfage_b;a<=slage_b;a++) //can still use age ranges for each survey
    {
      sigma2_Ip_md(y,a)=((1.-est_Ip_md(y,a))*est_Ip_md(y,a)+0.1/double(slage_b-sfage_b+1))/ESS_I_md;
    }//close age loop
  }//close year loop
  //cout << est_I_md << endl;
  //exit(1);
  mdssn=log(est_I_md);
  //ChesMMAP
  for(y=fmyear;y<=lmyear;y++)
  {
    for(t=1;t<=tstep;t++)
    {
      est_I_age_cm(t,y)=q_cm(t)*elem_prod(Nbay(t,y),ssel_cm(t));//year,timstep1, region1
     }//close tstep loop
  }//closing year loop
  //calculate total survey catch at age for each year
  for(t=1;t<=tstep;t++)
  {
    est_I_cm(t)=rowsum(est_I_age_cm(t));
  }//close tstep loop
  for(y=fmyear;y<=lmyear;y++)
  {//calculate proportions at age in the survey catch
    for(t=1;t<=tstep;t++)
    {
      est_Ip_cm(t,y)=est_I_age_cm(t,y)/est_I_cm(t,y);
      for(a=sfage_a;a<=slage_a;a++) //can still use age ranges for each survey
      {
        sigma2_Ip_cm(t,y,a)=((1.-est_Ip_cm(t,y,a))*est_Ip_cm(t,y,a)+0.1/double(slage_a-sfage_a+1))/ESS_I_cm;
      }//close age loop
    }//close tstep loop
  }//close year loop
  //cout << "chesmmap" << endl;
  chesmmap1=log(est_I_cm(1)); //sdreport ts 1
  chesmmap2=log(est_I_cm(2)); //sdreport ts2
  ////////////////////////
  //   Atlantic Coast   //
  ////////////////////////
  //CT List
  for(y=fmyear;y<=lmyear;y++)
  {
    for(t=1;t<=tstep;t++)
    {
      est_I_age_ct(t,y)=q_ct(t)*elem_prod(Ncoast(t,y),ssel_ct(t));//
      est_I_age_ct_err(t,y)=mod_age_err_a*est_I_age_ct(t,y); //catch at age with accounting for aging error
    }//close tstep loop
  }//closing year loop
  //calculate total survey catch at age for each year
  for(t=1;t<=tstep;t++)
  {
    est_I_ct(t)=rowsum(est_I_age_ct_err(t));//here s means survey, if we are hard coding this then that should disappear
  }//close tstep loop
  for(y=fmyear;y<=lmyear;y++)
  {//calculate proportions at age in the survey catch
    for(t=1;t<=tstep;t++)
    {
      est_Ip_ct(t,y)=est_I_age_ct_err(t,y)/est_I_ct(t,y);
      for(a=sfage_a;a<=slage_a;a++) //can still use age ranges for each survey
      {
        sigma2_Ip_ct(t,y,a)=((1.-est_Ip_ct(t,y,a))*est_Ip_ct(t,y,a)+0.1/double(slage_a-sfage_a+1))/ESS_I_ct;
      }//close age loop
    }//close tstep loop
  }//close year loop
  //cout << "ct" << endl;
  ctlist1=log(est_I_ct(1)); //sdreport ts 1
  ctlist2=log(est_I_ct(2)); //sdreport ts 1
  //NY Ocean Haul Survey
  for(y=fmyear;y<=lmyear;y++)
  {
    est_I_age_ny(y)=q_ny*elem_prod(Ncoast(2,y)(sfage_c,slage_c),ssel_ny);//
    est_I_age_ny_err(y)=mod_age_err_c*est_I_age_ny(y);
  }//closing year loop
  //calculate total survey catch at age for each year
  est_I_ny=rowsum(est_I_age_ny_err);//here s means survey, if we are hard coding this then that should disappear
  for(y=fmyear;y<=lmyear;y++)
  {//calculate proportions at age in the survey catch
    est_Ip_ny(y)=est_I_age_ny_err(y)/est_I_ny(y);
    for(a=sfage_c;a<=slage_c;a++) //can still use age ranges for each survey
    {
      sigma2_Ip_ny(y,a)=((1.-est_Ip_ny(y,a))*est_Ip_ny(y,a)+0.1/double(slage_c-sfage_c+1))/ESS_I_ny;
    }//close age loop
  }//close year loop
  //cout << "ny" << endl;
  nyohs=log(est_I_ny); //sdreport
  //NJ Bottom Trawl
  for(y=fmyear;y<=lmyear;y++)
  {
    est_I_age_nj(y)=q_nj*elem_prod(Ncoast(1,y)(sfage_b,slage_b),ssel_nj);//
    est_I_age_nj_err(y)=mod_age_err_b*est_I_age_nj(y);
  }//closing year loop
  //calculate total survey catch at age for each year
  est_I_nj=rowsum(est_I_age_nj_err);//here s means survey, if we are hard coding this then that should disappear
  for(y=fmyear;y<=lmyear;y++)
  {//calculate proportions at age in the survey catch
    est_Ip_nj(y)=est_I_age_nj_err(y)/est_I_nj(y);
    for(a=sfage_b;a<=slage_b;a++) //can still use age ranges for each survey
    {
      sigma2_Ip_nj(y,a)=((1.-est_Ip_nj(y,a))*est_Ip_nj(y,a)+0.1/double(slage_b-sfage_b+1))/ESS_I_nj;
    }//close age loop
  }//close year loop
  //cout << "nj" << endl;
  njbt=log(est_I_nj); //sdreport
  //DE SSN
  for(y=fmyear;y<=lmyear;y++)
  {
    //est_I_age_des(y)=q_des*elem_prod(Ncoast(1,y)(sfage_c,slage_c),ssel_des);//
    est_I_age_des(y)=q_des*elem_prod(N(2,2,1,y)(sfage_c,slage_c),ssel_des);//N(s,r,t,y,a)
    est_I_age_des_err(y)=mod_age_err_c*est_I_age_des(y);
  }//closing year loop
  //calculate total survey catch at age for each year
  est_I_des=rowsum(est_I_age_des_err);//here s means survey, if we are hard coding this then that should disappear
  for(y=fmyear;y<=lmyear;y++)
  {//calculate proportions at age in the survey catch
    est_Ip_des(y)=est_I_age_des_err(y)/est_I_des(y);
    for(a=sfage_c;a<=slage_c;a++) //can still use age ranges for each survey
    {
      sigma2_Ip_des(y,a)=((1.-est_Ip_des(y,a))*est_Ip_des(y,a)+0.1/double(slage_c-sfage_c+1))/ESS_I_des;
    }//close age loop
  }//close year loop
  //cout << "de ssn" << endl;
  dessn=log(est_I_des); //sdreport
  //DE 30
  for(y=fmyear;y<=lmyear;y++)
  {
    est_I_age_de30(y)=q_de30*elem_prod(Ncoast(2,y),ssel_de30);//
  }//closing year loop
  //calculate total survey catch at age for each year
  est_I_de30=rowsum(est_I_age_de30);//here s means survey, if we are hard coding this then that should disappear
  for(y=fmyear;y<=lmyear;y++)
  {//calculate proportions at age in the survey catch
    est_Ip_de30(y)=est_I_age_de30(y)/est_I_de30(y);
    for(a=sfage_a;a<=slage_a;a++) //can still use age ranges for each survey
    {
      sigma2_Ip_de30(y,a)=((1.-est_Ip_de30(y,a))*est_Ip_de30(y,a)+0.1/double(slage_a-sfage_a+1))/ESS_I_de30;
    }//close age loop
  }//close year loop
  de30=log(est_I_de30); //sdreport
  //cout << "de30" << endl;
  ///*********************************************
  ///            YOY and age 1 survey
  ///********************************************
  //age 1 surveys
  //cout << q_age1n << endl;
  for(y=fmyear;y<=lmyear;y++)
  {
    //NY age 1
    //for(t=1;t<=tstep;t++)
    //{
    for(o=1;o<=age1surv;o++)
    {
      t=2;
      //cout << "t=" << t << endl;
      est_I_age1n(o,y)=q_age1n*Ncoast(t,y,1);
      //cout << est_I_age1n(o,y) << endl;
      }//close o loop
    //}//close tstep loop
    //est_I_age1n(y)=q_age1n*Ncoast(2,y,1);//tstep=2, year, age=1
   // MD AGE 1
   t=1;
    est_I_age1m(y)=q_age1m*Nbay(t,y,1); //tstep=1, year, age=1
  }// close year looop
  //cout << "end age 1" << endl;
  //yoy surveys
  t=2;
  for(y=fmyear;y<=lmyear-1;y++)
  {
    //Cost survesy, NJYOY and NY YOY
    for(z=1;z<=yoysurv_coast;z++)  //second part of loop was missing "<" *MW
    {
      est_I_yoy_coast(z,y)=q_yoy_coast(z)*Ncoast(t,y+1,1);// time step 2?, year, age=1
     //cout << y << " " <<  z << " " << q_yoy_coast(z) << " " << Ncoast(2,y+1,1) << " " << est_I_yoy_coast(z,y) << endl;
    }//close z loop
    for(z=1;z<=yoysurv_bay;z++)
    {
       est_I_yoy_bay(z,y)=q_yoy_bay(z)*Nbay(t,y+1,1);//time step 2, year, age=1?
      //cout << y << " " <<  z << " " << q_yoy_coast(z) << " " << Ncoast(2,y+1,1) << " " << est_I_yoy_coast(z,y) << endl;
    }//close z loop
  }//close year loop
  //cout <<
  //exit(1);
  //cout << "finish yoy" << endl;
  //cout << "est I YOY coast " << endl << est_I_yoy_coast << endl << endl << endl;
  //cout << "est I YOY bay " << endl << est_I_yoy_bay << endl << endl << endl;
  //exit(1);
  oceanage1=log(est_I_age1n(1)); //sdreport
  bayage1=log(est_I_age1m); //sdreport
  for(y=fmyear;y<=lmyear-1;y++)
  {
    cbyoy(y)=log(est_I_yoy_bay(1,y)); //sdreport
    njyoy(y)=log(est_I_yoy_coast(1,y)); //sdreport
    nyyoy(y)=log(est_I_yoy_coast(2,y)); //sdreport
  }
}

void model_parameters::evaluate_likelihood(void)
{
  //lognormal likelihood for total catch
  for(t=1;t<=tstep;t++)
  {
    //cout << "t=" << t << " " << tstep << endl;
    for(r=1;r<=region;r++)
    {
      //cout << "r=" << r << endl;
      Lcatch(r,t)=lognorm_negLL(obs_C(r,t),est_region_C(r,t),C_var(r,t),fmyear,lmyear);
      //cout << est_region_C(r,t) << endl;
    }//close region loop
  }//close time step loop
   //cout << Lcatch << endl;
  //cout << "Lcatch" << endl;
  //exit(1);
  for(t=1;t<=tstep;t++)
  {  
    for(r=1;r<=region;r++)
    {
      Lcatchagecomp(r,t)=multinom_negLL(obs_Cp(r,t),est_region_Cp(r,t),sigma2_Cp(r,t),fage,lage,fmyear,lmyear);//multinomial for proportions
    }//close region loop
  }//close timestep loop
  //cout << "Lcatchagecomp" << endl;
  //lognormal likelkihood for indices of abundance
  //multinomial for age composition
  //Maryland
  Lindex_md=lognorm_negLL(obs_I_md,est_I_md,I_var_md,fmyear,lmyear);
  Lindexagecomp_md=multinom_negLL(obs_Ip_md,est_Ip_md,sigma2_Ip_md,sfage_b,slage_b,fmyear,lmyear);//change AGE
  //Lindex_md=0.;
  //Lindexagecomp_md=0.;
  //NY OHS
  Lindex_ny=lognorm_negLL(obs_I_ny,est_I_ny,I_var_ny,fmyear,lmyear);
  Lindexagecomp_ny=multinom_negLL(obs_Ip_ny,est_Ip_ny,sigma2_Ip_ny,sfage_c,slage_c,fmyear,lmyear);//change AGE
  //Lindex_ny=0.;
  //Lindexagecomp_ny=0.;
  //Nj BT
  Lindex_nj=lognorm_negLL(obs_I_nj,est_I_nj,I_var_nj,fmyear,lmyear);
  Lindexagecomp_nj=multinom_negLL(obs_Ip_nj,est_Ip_nj,sigma2_Ip_nj,sfage_b,slage_b,fmyear,lmyear);//change AGE
  //DE SSN
  Lindex_des=lognorm_negLL(obs_I_des,est_I_des,I_var_des,fmyear,lmyear);
  Lindexagecomp_des=multinom_negLL(obs_Ip_des,est_Ip_des,sigma2_Ip_des,sfage_c,slage_c,fmyear,lmyear);//change AGE
  //Lindex_des=0.;
  //Lindexagecomp_des=0.;
  //DE 30
  Lindex_de30=lognorm_negLL(obs_I_de30,est_I_de30,I_var_de30,fmyear,lmyear);
  Lindexagecomp_de30=multinom_negLL(obs_Ip_de30,est_Ip_de30,sigma2_Ip_de30,sfage_a,slage_a,fmyear,lmyear);//change AGE
  //cout << "obs i md"  << endl << obs_I_md << endl << "est i md" << endl << est_I_md << endl << "var" << endl << I_var_md << endl;
  //cout << "end one tstep surv" << endl;
  //cout << est_I_de30 << endl;
  for(t=1;t<=tstep;t++)
  {
    //CT LIST
    Lindex_ct(t)=lognorm_negLL(obs_I_ct(t),est_I_ct(t),I_var_ct(t),fmyear,lmyear);
    Lindexagecomp_ct(t)=multinom_negLL(obs_Ip_ct(t),est_Ip_ct(t),sigma2_Ip_ct(t),sfage_a,slage_a,fmyear,lmyear);//change AGEEE
    //chesMMAP
    Lindex_cm(t)=lognorm_negLL(obs_I_cm(t),est_I_cm(t),I_var_cm(t),fmyear,lmyear);
    Lindexagecomp_cm(t)=multinom_negLL(obs_Ip_cm(t),est_Ip_cm(t),sigma2_Ip_cm(t),sfage_a,slage_a,fmyear,lmyear);
  }
  //Lindex_ct(1)=0.;
  //Lindex_ct(2)=0.;
  //Lindexagecomp_ct(1)=0.;
  //Lindexagecomp_ct(2)=0.;
  /*
  cout << "ct " << endl;
  cout << "Obs IOA" << endl << obs_I_ct << endl;
  cout << "Est IOA" << endl << est_I_ct << endl;
  cout << "OBS CAA" << endl << obs_Ip_ct << endl;
  cout << "EST CAA" << endl << est_Ip_ct << endl;
  cout << "Lindex" << endl << Lindex_ct <<endl;
  cout << "Lindexagecomp" << endl << Lindexagecomp_ct << endl;
  exit(1);
  */
  //**************************
  //age 1 surveys
  //**************************
  //MD age 1
  Lage1index_md=lognorm_negLL(obs_I_age1m,est_I_age1m,I_var_age1m,fmyear,lmyear);// double check params
  //Lage1index_md=0.;
  //NY Age 1
  for(o=1;o<=age1surv;o++)
  {
    Lage1index_ny(o)=lognorm_negLL(obs_I_age1n(o),est_I_age1n(o),I_var_age1n(o),fmyear,lmyear);
  }
  //Lage1index_ny=0.;
  //cout << "end age 1 lik" << endl;
  //*********************
  //YOY surveys
  //********************
  for(z=1;z<=yoysurv_coast;z++)
  {
    Lyoyindex_coast(z)=lognormyoy_negLL(obs_I_yoy_coast(z),est_I_yoy_coast(z),I_var_yoy_coast(z),fmyear,lmyear-1);//lmyear-1 because cannot caluclate yoy in the last year
    //Lyoyindex_coast(z)=0.;
  }
  for(z=1;z<=yoysurv_bay;z++)
  {
    Lyoyindex_bay(z)=lognormyoy_negLL(obs_I_yoy_bay(z),est_I_yoy_bay(z),I_var_yoy_bay(z),fmyear,lmyear-1);
    //Lyoyindex_bay(z)=0.;
  }
  //cout << "end yoy like" << endl;
  //cout << Lyoyindex_bay << endl;
  //exit(1);
  //pen_N0_dev=0.5*norm2(log_N0_devs/0.49);
  //pen_f2sel=0.5*norm2((log_fs(2,1)-mean(log_fs(2,1))))/0.25;
  if(!last_phase())
  {
    pen_F=0.5*square(sum(log_F)-log(0.3));//assuming normal penalty for log(F)
  }
  //pen_r=
  pen_rdev=norm2(log_Rdevs1)+norm2(log_Rdevs2);//recruitment deviation penalty; assumes logscale variable of mean recruiment has SD 1 and mean of 0
  pen_fdev=0.;
  for(y=fmyear+1;y<=lmyear;y++)
  {
    pen_fdev+=1./square(0.5)*square(log_Fdevs_r1t1(y)-log_Fdevs_r1t1(y-1)); //penalizing the model for differences in F from one year to the next year
    pen_fdev+=1./square(0.5)*square(log_Fdevs_r1t2(y)-log_Fdevs_r1t2(y-1));
    pen_fdev+=1./square(0.5)*square(log_Fdevs_r2t1(y)-log_Fdevs_r2t1(y-1));
    pen_fdev+=1./square(0.5)*square(log_Fdevs_r2t2(y)-log_Fdevs_r2t2(y-1));
    //square 0.5 suggests that the penalty should have a SD of 0.5 on a normal distribution
  }
  //penalty for occupancy probabilities deviating from priors
  pen_prop=0.;//setting penalty to 0 to start
  for(ts=1;ts<=tstep;ts++)
  {
    for(a=fage+1;a<=lage;a++)
    {
      pen_prop+=0.5*square((log_prop_bay(ts,a)-prop_bay(ts,a))/log_sd_bay(ts,a));//adding for bay fish in the bay
       // use_aco_prop is the  switch to multiply by 0 if dat file indicates
      //pen_prop+=0.5*square((log(acoustic_prop(ts,a)+0.01)-prop_bay(ts,a))/log(acoustic_prop_sd(ts,a)+0.01))*double(use_aco_prop);//add penalty for the acousitc occupancy probabilities. use_aco_prop is a switch to tell if this should be included in the likelihood or not
      //pen_prop+=beta_mig(prop(1,ts,2),alpha_bay(ts),beta_bay(ts)); //feeding in the  occupancy probabilites, alpha, and beta parameters for stock 1 (Bay) in the Coast, for each time step
    }//close age loop
  }//close ts loop
  //cout << pen_prop << endl;
  //exit(1);
  //add penalties for atl coast stock to pen_prop
  for(ts=1;ts<=tstep;ts++)
  {
    for(a=fage+1;a<=lage;a++)
    {
      pen_prop+=0.5*square((log_prop_coast(ts,a)-prop_coast(ts,a))/log_sd_coast(ts,a));//adding for coast fish in the coast
    }//close age loop
  }//close ts loop
  //add acoustic occupancy proababilities to pen_prop
  pen_prop_aco=0.0;
  for(ts=1;ts<=tstep;ts++)
  {
    for(a=fage+1;a<=lage;a++)
    {
        //if(use_aco_prop==2)//if region = chesapeake bay
      //{
       // pen_prop+=0.0;
      //}
      //else
      //{
        if(acoustic_prop_sd(ts,a)!=-99)
        {
          pen_prop_aco+=0.5*square((log((1-acoustic_prop(ts,a))+0.01)-prop_bay(ts,a))/acoustic_prop_sd(ts,a)); //prop sd is entered on the log scale already, techinally the CV (log(SD/mean+0.01)
        }
      //}//close else
    }
  }
  //cout << "with added penalty" << endl << pen_prop << endl;
  //cout << "log prop coast" << endl << log_prop_coast << endl;
  //cout << "prop coast" << endl << prop_coast << endl;
  //exit(1);
  //penalty for proportion of CB fish in the coast
  //pen_prop_bay=-765.*0.65*log(0.01+sum(N(1,2,2,lmyear)(4,lage))/sum(N(1,2,2,lmyear)(4,lage)+N(2,2,2,lmyear)(4,lage)));//binomial log likelihood, taking neg for NLL
  //765 was ESS
  //***removed first two from Hasegawa et al. 2022 because outside of model time frame
  //pen_prop_bay=-1000.*0.65*log(0.01+sum(est_C_age_err(1,2,2,lmyear)(4,lage))/sum(est_C_age_err(1,2,2,lmyear)(4,lage)+est_C_age_err(2,2,2,lmyear)(4,lage)));//binomial log likelihood, taking neg for NLL
  //pen_prop_bay+=-1000.*0.35*log(0.01+sum(est_C_age_err(2,2,2,lmyear)(4,lage))/sum(est_C_age_err(1,2,2,lmyear)(4,lage)+est_C_age_err(2,2,2,lmyear)(4,lage)));//binomial log likelihood, taking neg for NLL
  //next two penalties are from Kneebone et al. 2012
         //tagged 159 fish, no information on how many tripes
         //changed ESS from 1000 to 159, 3/25/24
  pen_prop_bay=-159.*0.61*log(0.01+sum(est_C_age_err(1,2,2,2010)(4,lage))/sum(est_C_age_err(1,2,2,2010)(4,lage)+est_C_age_err(2,2,2,2010)(4,lage)));//binomial log likelihood, taking neg for NLL
  pen_prop_bay+=-159.*0.39*log(0.01+sum(est_C_age_err(2,2,2,2010)(4,lage))/sum(est_C_age_err(1,2,2,2010)(4,lage)+est_C_age_err(2,2,2,2010)(4,lage)));
  //  765= ESS, 0.5 = observed proportion of chesapeake bay fish in the coast during time step 2; estimated proporortion
  //0.01, added constant so it's not the log of 0
  //pen_prop_bay=0.0;
  //cout << "pen prop" << endl << pen_prop_bay << endl;
  //cout << "est c age err" << endl << est_C_age << endl;
  //exit(1);
  //cout << log_sf3_cb(1,2) << endl;
  //exit(1);
  //prior for descending limb of chesapeake bay fsel
   pen_cb_sel=//0.5*square((log_sf1_cb(1,2)-0.)/0.5)+ //ascending param slope for the first timeblock
             0.5*square((log_sf2_cb(1,2)-1.098)/0.5)+ //ascending param 50% sel
             0.5*square((log_sf3_cb(1,1)-0.)/0.25)+ //descending param slope for the first timeblock
             0.5*square((log_sf4_cb(1,1)-4.)/0.25)+//+//descending param 50% selc for the first timeblock
             0.5*square((log_sf3_cb(1,2)-0.)/0.25)+ //descending param slope for the first timeblock
             0.5*square((log_sf4_cb(1,2)-4.)/0.25);//+//descending param 50% selc for the first timeblock
             //0.5*square((log_sf3_cb(2,2)-0.)/0.5)+//descending param slope for the 2nd timeblock
             //0.5*square((log_sf4_cb(2,2)-4.)/0.5)+//descending param 50% selc for the 2nd timeblock
             //0.5*square((log_sf3_cb(3,2)-0.)/0.5)+//descending param slope for the 3rd timeblock
             //0.5*square((log_sf4_cb(3,2)-4.)/0.5)+ //descending params for fsel in the 3rdd timeblock
             //0.5*square((log_sf3_cb(4,2)-0.)/0.5)+//descending param slope for the 4th timeblock
             //0.5*square((log_sf4_cb(4,2)-4.)/0.5); //descending params for fsel in the 4th timeblock
  //pen_cb_sel=0.;
  pen_sf=0.0;
  if(log_sf4_cb(3,2)<log_sf2_cb(3,2))
  {
    pen_sf=10.*square(log_sf4_cb(3,2)-log_sf2_cb(3,2));
  }
  pen_sf_ct=0.0;
  if(log_ssf4_ct(1)<log_ssf2_ct(1))
  {
    pen_sf_ct=10.*square(log_ssf4_ct(1)-log_ssf2_ct(1));
  }
  if(log_ssf4_ct(2)<log_ssf2_ct(2))
  {
    pen_sf_ct+=10.*square(log_ssf4_ct(2)-log_ssf2_ct(2));
  }
  pen_feq=0.5*square((mfexp(log_Feq(1))-0.3)/0.1);
  //add all the components of the negative log likelihood together
  //neg_LL=sum(Lcatch)+sum(Lcatchagecomp)+sum(Lindex_a)+sum(Lindex_b)+sum(Lindex_c)+sum(Lindexagecomp_a)+sum(Lindexagecomp_b)+sum(Lindexagecomp_c)+sum(Lage1index)+sum(Lyoyindex)+Lfsel+Lssel_a+Lssel_b+Lssel_c+pen_N0_dev+pen_f2sel+pen_F; 
  //cout << "Lindex nj" << " " << Lindex_nj << endl;
  neg_LL=sum(Lcatch)+sum(Lcatchagecomp)+Lindex_md+sum(Lindex_cm)+Lindex_ny+Lindex_nj+sum(Lindex_ct)+Lindex_des+
          Lindex_de30+Lindexagecomp_md+sum(Lindexagecomp_cm)+Lindexagecomp_ny+Lindexagecomp_nj+sum(Lindexagecomp_ct)+
          Lindexagecomp_des+Lindexagecomp_de30+Lage1index_md+sum(Lage1index_ny)+sum(Lyoyindex_bay)+sum(Lyoyindex_coast)+
          pen_F+pen_prop+pen_prop_aco+pen_prop_bay+pen_cb_sel+pen_rdev+pen_fdev+pen_feq+pen_sf+pen_sf_ct;//sum(Lfsel)+Lssel_md+sum(Lssel_cm)+Lssel_ny+Lssel_nj+sum(Lssel_ct)+Lssel_des+Lssel_de30
  /*
  cout << "neg_LL " << endl << neg_LL << " = " << sum(Lcatch) << " sum(Lcatch)" << " + " << sum(Lcatchagecomp) << " sum (Lcatchagecomp" <<  " + " << Lindex_md << "Lindex_md" << " + " <<
  sum(Lindex_cm) << endl << " L index cm + " << Lindex_ny << " L Index ny + " << Lindex_nj << " Lindex nj +" << sum(Lindex_ct) << " L index ct + " << endl << Lindex_des << " L index des + "
  << Lindex_de30<< " L Index de30 + " << Lindexagecomp_md << " L index agecomp md + " << sum(Lindexagecomp_cm) << "L index agecomp cm + " << Lindexagecomp_ny << endl << "Lindex agecomp ny + "
  << Lindexagecomp_nj << "L index agecomp nj + " << sum(Lindexagecomp_ct) << "L index agecomp ct + " << Lindexagecomp_des << "L index agecomp des + " << Lindexagecomp_de30 << "L index agecomp de30 + "
  << Lage1index_md << "Lindex agecomp md + " << sum(Lage1index_ny) << "L age1 index ny + " << Lyoyindex_bay << "Lyoyindex bay + " << endl << sum(Lyoyindex_coast) << "Lyoyindex coast + " << sum(Lfsel) <<
  "Lfsel + " << Lssel_md << "L ssel md + " << sum(Lssel_cm) << "Lssel cm + " << Lssel_ny << " Lssel ny + " << Lssel_nj << "Lssel nj + " << sum(Lssel_ct) << "lsselct + " << Lssel_des << " Lssel des + " <<
  Lssel_de30 << "Lssel de30 + " << pen_F << " pen f + " << pen_prop << " pen prop + " << pen_prop_bay << "pen prop bay + " << pen_cb_sel << " pen cb sel" << pen_rdev << " + pen rdev " 
  << pen_fdev << " + pen f dev " << endl <<endl;
  //exit(1);
  //cout << "end neg_ll" << endl;
  //exit(1);
  */
}

dvar_vector model_parameters::migrate(dvar_vector N,dvar_vector P)
{
  return(elem_prod(N,P));
}

dvariable model_parameters::beta_mig(dvar_vector occ_prob, dvector alpha, dvector beta)
{
  dvariable betaLL;
  betaLL=0.0;
  for(a=fage+1;a<=lage;a++)
  {
    betaLL+=-(alpha(a)-1)*log(occ_prob(a))-(beta(a)-1)*log(1-occ_prob(a));
  }
  return(betaLL);
}

dvariable model_parameters::lognorm_negLL(dvector obsI, dvar_vector estI, dvector Ivar, data_int fmyear, data_int lmyear)
{
  dvariable negLL;
  negLL=0.0;
  for(y=fmyear;y<=lmyear;y++)
  {
    if(obsI(y)!=-99)
    {
      negLL+=square((log(obsI(y)+0.01)-log(estI(y)+0.01)))/(2.*Ivar(y));
    }
  }
  //negLL=norm2(elem_div(log(obsI(fmyear,lmyear))-log(estI),(2.*Ivar)));
  return(negLL);
}

dvariable model_parameters::lognormyoy_negLL(dvector obsI, dvar_vector estI, dvector Ivar, data_int fmyear, data_int lmyear)
{
  dvariable negLL;
  negLL=0.0;
  for(y=fmyear;y<=lmyear;y++)
  {
    if(obsI(y)!=-99)
    {
      negLL+=square((log(obsI(y))-log(estI(y))))/(2.*Ivar(y));
    }
  }
  //negLL=norm2(elem_div(log(obsI(fmyear,lmyear))-log(estI),(2.*Ivar)));
  return(negLL);
}

dvariable model_parameters::lognormyoy_negLL(dvector obsI, dvar_vector estI, dvector Ivar, data_int fmyear, int lmyear)
{
  dvariable negLL;
  negLL=0.0;
  for(y=fmyear;y<=lmyear;y++)
  {
    if(obsI(y)!=-99)
    {
      negLL+=square((log(obsI(y))-log(estI(y))))/(2.*Ivar(y));
    }
  }
  //negLL=norm2(elem_div(log(obsI(fmyear,lmyear))-log(estI),(2.*Ivar)));
  return(negLL);
}

dvariable model_parameters::multinom_negLL(data_matrix obsP, named_dvar_matrix estP, named_dvar_matrix sigma2, int fage, int lage, int fyear, int lyear)
{
  dvariable negLL; //declaring new variable called negLL
  negLL=0.0;
  for(y=fyear;y<=lyear;y++)
  {
    for(a=fage;a<=lage;a++)
    {
      //L4+=-ESS_C*obs_Ip(y)*log(est_Ip(y)+0.001);
      if(obsP(y,a)!=-99)//runs likelihood if value is not equal to -99, accounting for missing values
      {
        negLL+=-log(mfexp(.5*(-square(obsP(y,a)-estP(y,a))/sigma2(y,a)))+0.01);
      }
    }
  }
  return(negLL);
}

dvariable model_parameters::multinom_negLL(dmatrix obsP, dvar_matrix estP, dvar_matrix sigma2, int fage, int lage, int fyear, int lyear)
{
  dvariable negLL; //declaring new variable called negLL
  negLL=0.0;
  for(y=fyear;y<=lyear;y++)
  {
    for(a=fage;a<=lage;a++)
    {
      //L4+=-ESS_C*obs_Ip(y)*log(est_Ip(y)+0.001);
      if(obsP(y,a)!=-99)//runs likelihood if value is not equal to -99, accounting for missing values
      {
        negLL+=-log(mfexp(.5*(-square(obsP(y,a)-estP(y,a))/sigma2(y,a)))+0.01);
      }
    }
  }
  return(negLL);
}

dvariable model_parameters::multinom_negLL(dmatrix obsP, dvar_matrix estP, dvar_matrix sigma2, data_int fage, data_int lage, int fyear, int lyear)
{
  dvariable negLL; //declaring new variable called negLL
  negLL=0.0;
  for(y=fyear;y<=lyear;y++)
  {
    for(a=fage;a<=lage;a++)
    {
      //L4+=-ESS_C*obs_Ip(y)*log(est_Ip(y)+0.001);
      if(obsP(y,a)!=-99)//runs likelihood if value is not equal to -99, accounting for missing values
      {
        negLL+=-log(mfexp(.5*(-square(obsP(y,a)-estP(y,a))/sigma2(y,a)))+0.01);
      }
    }
  }
  return(negLL);
}

dvariable model_parameters::multinom_negLL(data_matrix obsP, named_dvar_matrix estP, named_dvar_matrix sigma2, data_int fage, data_int lage, int fyear, int lyear)
{
  dvariable negLL; //declaring new variable called negLL
  negLL=0.0;
  for(y=fyear;y<=lyear;y++)
  {
    for(a=fage;a<=lage;a++)
    {
      //L4+=-ESS_C*obs_Ip(y)*log(est_Ip(y)+0.001);
      if(obsP(y,a)!=-99)//runs likelihood if value is not equal to -99, accounting for missing values
      {
        negLL+=-log(mfexp(.5*(-square(obsP(y,a)-estP(y,a))/sigma2(y,a)))+0.01);
      }
    }
  }
  return(negLL);
  /*
}

dvariable model_parameters::logisticnorm_negLL(dvector obsP, dvar_vector estP, prevariable sigma2, data_int fage, data_int lage, data_int fyear, data_int lyear, double ESS,  prevariable phi, int B )
{
  dvariable negLL; //declaring variable called negLL
  //ivector B(fyear,lyear); //max age with observerd proportion > 0
  dvar_matrix V(fage,B-1,fage,B-1); //declaring variable called V, defined as Bins*variance^(2*(bins-1))
  dvar_vector w(fage,B-1); //declaring variable called w, defined as a vector with length bins-1
  //dvar_vector W(fyear,lyear); //declaring capital W, defined as squre root fo meanNy/Ny
  dvar_matrix Vinv(fage,B-1,fage,B-1);//calculating B inverse
  //dvar_vector (1,B);
  //double b;
  //b=double(B);
  int a1;
  int a2;
 // W=sqrt(ESS/double(lyear-fyear+1));
  negLL=0.0; //starting calculations for logistic normal likelihood
  for(a1=fage;a1<=B-1;a1++)
    {
      V(a1,a1)=sigma2; //the diagonals, equation from Francis 2014, pow(base,root) is raised to a power
      for(a2=fage;a2<=B-1;a2++)
      {
        if(a2!=a1)
        {
          V(a1,a2)=pow(phi,abs(double(a1-a2)))*V(a1,a1); // correlation * variance
        }
       //cout << a1 << " " << a2 << endl;
      }
    }
    Vinv=inv(V);
    for(a=fage;a<=B-1;a++)
        {
         if(obsP(a)!=-99)
         {
           if(obsP(a)>0 && estP(a)>0 && estP(B)>0)//|| is 'or'
           {
            w(a)=log(obsP(a)/obsP(B))-log(estP(a)/estP(B)); //equatioon from Fisch et al 2021, and francis 2014
           }//close if statement
           else
           {
             w(a)=0.0;
           }
         }
    }//close age loop
    //cout << "w " << w << endl;
    //cout << "Obs Pa " << obsP << endl;
    //cout << "Est Pa " << estP << endl;
    for(a=fage;a<=B-1;a++)
    {
    if(obsP(a)>0 && estP(a)>0 && estP(B)>0)
      //if(obsP(a)>0)
      {
        negLL+=0.5*double(B-1)*log(2.0*PI)+log(obsP(a))+0.5*log(det(V)); //equation from Francis 2014 and Fisch 2021
      }//close if statement
      //cout << "negLL " << negLL << endl;
    }//close age loop
    negLL+=0.5*w*(Vinv)*w;//last term of likelihood equation
  return(negLL);
  */
}

void model_parameters::report(const dvector& gradients)
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
  //The report section is used to write output to the standard output "filename.rep"
  report << "observed catch" << endl << obs_C << endl;
  report << "estimated catch" << endl << est_region_C << endl;
  //report << "observed index" << endl << obs_I << endl;
  //report << "estimated index" << endl << est_I << endl;
  report << "fishery selectivity" << endl << fsel << endl;
  //report << "survey selectivity" << endl << ssel << endl;
  report << "abundance" << endl << N << endl;
  report << "fishing mortality" << endl << F << endl;
  report << "biomass" << endl << B << endl;
  report << "SSB" << endl << SSB << endl;
  ofstream catout("catch.txt");
  {
    catout << "year timestep region logobs logpred " << endl;
    for(y=fmyear;y<=lmyear;y++)
    {
      for(t=1;t<=tstep;t++)
      {
       for(r=1;r<=region;r++)
       {
         catout << y << " " << t << " " << r << " "  << log(obs_C(r,t,y)+0.00001) << " " << log(est_region_C(r,t,y)+0.00001) << endl;
        }//close region
      }//close tstep
    }//close yr
   }//close ofstream
    ofstream ioaout("ioa.txt");
    {
        ioaout << "agegroup survey year timestep obsioa predioa sd" << endl;
        for(y=fmyear;y<=lmyear;y++)
        {
          ioaout << "1-15" << " " <<  "DE30" << " " << y << " " << "2" << " "<< log(obs_I_de30(y)+0.00001) << " " << log(est_I_de30(y)+0.00001) << " " << obs_I_CV_de30(y) <<  endl;
          ioaout << "2-15" << " " << "MD" << " " << y << " " << "1" << " "<< log(obs_I_md(y)+0.00001) << " " << log(est_I_md(y)+0.00001) << " " << obs_I_CV_md(y) << endl;
          ioaout << "2-15" << " " << "NJ" << " " << y << " " << "1" << " "<< log(obs_I_nj(y)+0.00001) << " " << log(est_I_nj(y)+0.00001) << " " << obs_I_CV_nj(y) << endl;
          ioaout << "2-13" << " " << "NY" << " " << y << " " << "2" << " "<< log(obs_I_ny(y)+0.00001) << " " << log(est_I_ny(y)+0.00001) << " " << obs_I_CV_ny(y) << endl;
          ioaout << "2-13" << " " << "DESSN" << " " << y << " " << "1" << " "<< log(obs_I_des(y)+0.00001) << " " << log(est_I_des(y)+0.00001) << " " << obs_I_CV_des(y) << endl;
          for(t=1;t<=tstep;t++)
          {
             ioaout << "1-15" << " " << "cm" << " " << y << " " << t << " "<< log(obs_I_cm(t,y)+0.00001) << " " << log(est_I_cm(t,y)+0.00001) << " " << obs_I_CV_cm(t,y) << endl;
            ioaout << "1-15" <<  " " << "ct" << " " << y << " " << t << " "<< log(obs_I_ct(t,y)+0.00001) << " " << log(est_I_ct(t,y)+0.00001) << " " << obs_I_CV_ct(t,y) << endl;
          }//close timestep
       }//close year
     }//close ofstream
   ofstream obscaa("ocaa.txt");//obs fishery catch at age
   {
     obscaa << " region timestep year age obscat obscaa estcaa standardresid" << endl;
     for(y=fmyear;y<=lmyear;y++)
     {
       for(t=1;t<=tstep;t++)
       {
           for(a=fage;a<=lage;a++)
           {
             for(r=1;r<=region;r++)
             { 
               obscaa << " " << r << " " << t << " " << y << " " << a << " " << obs_C(r,t,y)*obs_Cp(r,t,y,a) << " " <<  obs_Cp(r,t,y,a) << " " << est_region_Cp(r,t,y,a) << " " <<  log((obs_Cp(r,t,y,a)-est_region_Cp(r,t,y,a))+0.001)/sqrt(log(obs_Cp(r,t,y,a)+0.001)) << endl;
             }//close region
           }//close age      
     }//close tstep
     }//close yr
   }//close ofstream
   ofstream obssaa("osaa.txt");
   {
     obssaa << "agegroup survey year timestep age obscat obssaa estsaa standardresid" << endl;
     for(y=fmyear;y<=lmyear;y++)
     {
       for(a=sfage_a;a<=slage_a;a++)
       {
         if(obs_Ip_de30(y,a)!=-99) obssaa << "1-15 " << " " << "de30"  << " " << y << " " << "2" << " " << a << " " << obs_I_de30(y)*obs_Ip_de30(y,a) << " " <<  obs_Ip_de30(y,a) << " " <<  est_Ip_de30(y,a) << " " << (obs_Ip_de30(y,a)-est_Ip_de30(y,a))/sqrt(I_var_de30(y)) << endl;
       }//close age
     }//close year
     //obssaa << "agegroup survey year timestep age obssaa estsaa standardresid" << endl;
     for(y=fmyear;y<=lmyear;y++)
     {
       for(a=sfage_a;a<=slage_a;a++)
       {
       for(t=1;t<=tstep;t++)
         {
           if(obs_Ip_cm(t,y,a)!=-99) obssaa << "1-15" << " " << "cm"  << " " << y << " " << t << " " << a << " " << obs_I_cm(t,y)*obs_Ip_cm(t,y,a) << " " <<  obs_Ip_cm(t,y,a) << " " <<  est_Ip_cm(t,y,a) << " " << (obs_Ip_cm(t,y,a)-est_Ip_cm(t,y,a))/sqrt(I_var_cm(t,y)) << endl;
           if(obs_Ip_ct(t,y,a)!=-99) obssaa << "1-15" << " " << "ct"  << " " << y << " " << t << " " << a << " " << obs_I_ct(t,y)*obs_Ip_ct(t,y,a) << " " <<  obs_Ip_ct(t,y,a) << " " <<  est_Ip_ct(t,y,a) << " " << (obs_Ip_ct(t,y,a)-est_Ip_ct(t,y,a))/sqrt(I_var_ct(t,y)) << endl;
         }//close tstep
       }//close age
     }//close year
    for(y=fmyear;y<=lmyear;y++)
    {
      for(a=sfage_b;a<=slage_b;a++)
      {
         if(obs_Ip_nj(y,a)!=-99) obssaa << "2-15" << " " << "nj"  << " " << y <<  " " << "1" << " " << a << " " << obs_I_nj(y)*obs_Ip_nj(y,a) << " " << obs_Ip_nj(y,a) << " " <<  est_Ip_nj(y,a) << " " << (obs_Ip_nj(y,a)-est_Ip_nj(y,a))/sqrt(I_var_nj(y)) << endl;
         if(obs_Ip_md(y,a)!=-99) obssaa << "2-15" << " " << "md"  << " " << y << " " << "1" << " " << a << " " <<  obs_I_md(y)*obs_Ip_md(y,a) << " " << obs_Ip_md(y,a) << " " <<  est_Ip_md(y,a) << " " << (obs_Ip_md(y,a)-est_Ip_md(y,a))/sqrt(I_var_md(y)) << endl;
      }//close age
    }//close year
    for(y=fmyear;y<=lmyear;y++)
    {
      for(a=sfage_c;a<=slage_c;a++)
      {
         if(obs_Ip_ny(y,a)!=-99) obssaa << "2-13" << " " << "ny"  << " " << y  << " " << "2" << " " << a << " " << obs_I_ny(y)*obs_Ip_ny(y,a) << " " << obs_Ip_ny(y,a) << " " <<  est_Ip_ny(y,a) << " " << (obs_Ip_ny(y,a)-est_Ip_ny(y,a))/sqrt(I_var_ny(y)) << endl;
         if(obs_Ip_des(y,a)!=-99) obssaa << "2-13" << " " << "des"  << " " << y << " " << "1" << " " << a << " " << obs_I_des(y)*obs_Ip_des(y,a) << " " << obs_Ip_des(y,a) << " " <<  est_Ip_des(y,a) << " " << (obs_Ip_des(y,a)-est_Ip_des(y,a))/sqrt(I_var_des(y)) << endl;
       }//close age
     }//close year
   }//close of stream
   ofstream fselout("fsel.txt");//fsel at age
   {
     fselout << "region timestep year age fsel" << endl;
     for(r=1;r<=region;r++)
     {
       for(y=fmyear;y<=lmyear;y++)
       {
         for(t=1;t<=tstep;t++)
         {
           for(a=fage;a<=lage;a++)
           {
             fselout << r  << " " << t << " " << y << " " << a << " " << log(fsel(r,t,y,a)) << endl;
           }//age
         }//close tstep
       }//clse year
     }//closeregion
   }//close ofstream
   ofstream sselout("ssel.txt");//survey selectivty
   {
     sselout << "agegroup survey timestep age ssel" << endl;
     for(a=sfage_a;a<=slage_a;a++)
     {
       sselout << "1-15 " << "DE30" << " " << "2" << " "<< a << " " << log(ssel_de30(a)) << endl;
       for(t=1;t<=tstep;t++)
       {
         sselout << "1-15 " << "cm " << " " << t << " "<< a << " " << log(ssel_cm(t,a)) << endl;
         sselout << "1-15 " << "ct " << " " << t << " "<< a << " " << log(ssel_ct(t,a)) << endl;
         }//close tstep
     }//close a age
    for(a=sfage_b;a<=slage_b;a++)
    {
     sselout << "2-15 " << "NJ" << " " << "2" << " " << a << " " << log(ssel_nj(a)) << endl;
     sselout << "2-15 " << "MD"<< " " << "1" << " " << a << " " << log(ssel_md(a)) << endl;
    }
   for(a=sfage_c;a<=slage_c;a++)
   {
     sselout << "2-13 " << "NY" << " " << "2" << " " << a << " " << log(ssel_ny(a)) << endl;
     sselout << "2-13 " << "DESSN" << " " << "1" << " " << a << " " << log(ssel_des(a)) << endl;
    }
   }
   ofstream obsage1("age1.txt");
   {
     obsage1 << "survey year timestep obs est sd" << endl;
     for(y=fmyear;y<=lmyear;y++)
     {
       if(obs_I_age1m(y)!=-99) obsage1 << "md" << " "  << y << " " << "2" << " "  <<  log(obs_I_age1m(y)+0.001) << " " <<  log(est_I_age1m(y)+0.001) << " " << obs_I_age1m_CV(y) <<  endl;
       for(o=1;o<=age1surv;o++)
       {
           if(obs_I_age1n(o,y)!=-99) obsage1 << "ny"  << " " << y << " " <<  "2" << " "  <<  log(obs_I_age1n(o,y)+0.001) << " " <<  log(est_I_age1n(o,y)+0.001) << " " << obs_I_age1n_CV(o,y) << endl;
       }
     }
   }
   ofstream obsyoy("yoy.txt");
   {
     obsyoy << "region survey year obsyoy estyoy sd" << endl;
     for(y=fmyear;y<=lmyear;y++)
     {
      for(z=1;z<=yoysurv_coast;z++)
      {
        if(obs_I_yoy_coast(z,y)!=-99) obsyoy << "coast" << " " <<  z  << " " << y << " " <<  log(obs_I_yoy_coast(z,y)+0.001) << " " <<  log(est_I_yoy_coast(z,y)+0.001) << " " << obs_I_yoy_CV_coast(z,y) << endl;
       }//close z loop
      for(z=1;z<=yoysurv_bay;z++)
      {
        if(obs_I_yoy_bay(z,y)!=-99) obsyoy  << "bay" << " " << z << " " << y << " " <<  log(obs_I_yoy_bay(z,y)+0.001) << " " <<  log(est_I_yoy_bay(z,y)+0.001) << " " << obs_I_yoy_CV_bay(z,y) <<  endl;
   }//close z loop
     }//close year loop
   }//close of stream
   ofstream pop("pop.txt");// pop across years
   {
      pop << "stock region timstep year age  pop " << endl; 
      for(y=fmyear;y<=lmyear;y++)
      {
        for(t=1;t<=tstep;t++)
        {
          for(a=fage;a<=lage;a++)
          {
            for(r=1;r<=region;r++)
            {
              for(s=1;s<=stock;s++)
              {
               pop << " " << s << " " << r << " " <<  t << " " <<  y << " "  << a << " "  << N(s,r,t,y,a) << endl;
              }//close stock loop
            }//close region
         }//close age
      }//close tstep
     }//close year
   }
   ofstream mort("f.txt");// pop across years
   {
     mort << "region timestep year age  f " << endl; 
     for(r=1;r<=region;r++)
      {
           for(y=fmyear;y<=lmyear;y++)
           {
              for(t=1;t<=tstep;t++)
              {
                  for(a=fage;a<=lage;a++)
                  {
                  mort << " " << r << " " << t << " " << y  << " " << a << " " << F(r,t,y,a) << endl;
                  }//close age
              }//close tstep 
            }//close year
        }//close region loop
    }//close ofstream
   ofstream weightf("weightf.txt");
   {
     weightf << "region timestep year age  f " << endl;
          for(r=1;r<=region;r++)
      {
           for(y=fmyear;y<=lmyear;y++)
           {
              for(ts=1;ts<=tstep;ts++)
              {
                  for(a=7;a<=lage;a++)
                  {
                  weightf << " " << r << " " << ts << " " << y  << " " << a << " " << Freg_7plus(r,ts,y,a) << endl;
                  }//close age
              }//close tstep 
            }//close year
        }//close region loop
   }//close ofstream
   ofstream stockF("stockF.txt");
   {
     stockF << "year stock f" << endl;
     for(y=fmyear;y<=lmyear;y++)
     {
       for(s=1;s<=stock;s++)
       {
         stockF << " " << y << " " << s << " " << Fplus(s,y) <<  endl;
       }//close stock loop
     }//close year loop
   }//close ofstream
   ofstream like("lik.txt");
   {
     like << "name likelihood" << endl;
     like << " " << "sumLcatch_bay" << " " << Lcatch(1,1)+Lcatch(1,2) << endl;//adding for both timesteps
     like << " " << "sumLcatch_coast" << " " <<Lcatch(2,1)+Lcatch(2,2) << endl;
     like << " " << "sumLcatchagecomp_bay" << " " << Lcatchagecomp(1,1)+Lcatchagecomp(1,2) << endl;
     like << " " << "sumLcatchagecomp_coast" << " " << Lcatchagecomp(2,1)+Lcatchagecomp(2,2) << endl;   
     like << " " << "Lindexmd" << " " << Lindex_md << endl;
     like << " " << "Lindexcm " << " " << sum(Lindex_cm) << endl;
     like << " " << "Lindex_ny " << " " << Lindex_ny << endl;
     like << " " << "Lindex_nj" << " " << Lindex_nj << endl;
     like << " " << "Lindex_des" << " " << Lindex_des << endl;
     like << " " << "Lindex_de30" << " " << Lindex_de30 << endl;
     like << " " << "sum(Lindex_ct)" << " " << sum(Lindex_ct) << endl;
     like << " " << "Lindexagecomp_md" << " " << Lindexagecomp_md << endl;
     like << " " << "Lindexagecomp_cm" << " " << sum(Lindexagecomp_cm) << endl;
     like << " " << "Lindexagecomp_ny" << " " << Lindexagecomp_ny << endl;
     like << " " << "Lindexagecomp_nj" << " " << Lindexagecomp_nj << endl;
     like << " " << "Lindexagecomp_ct" << " " << sum(Lindexagecomp_ct) << endl;
     like << " " << "Lindexagecomp_des" << " " << Lindexagecomp_des << endl;
     like << " " << "Lindexagecomp_de30" << " " << Lindexagecomp_de30 << endl;
     like << " " << "Lindexage1md" << " " << Lage1index_md << endl;
     like << " " << "Lindexage1ny" << " " << Lage1index_ny << endl;
     like << " " << "Lindexcbyoy" << " " << sum(Lyoyindex_bay) << endl;
     like << " " << "Lindexacyoy" << " " << sum(Lyoyindex_coast) << endl;
     like << " " << "Pen_f" << " " << pen_F << endl;
     like << " " << "pen_prop_bay" << " " << pen_prop_bay << endl;
     like << " " << "pen_prop" << " " << pen_prop << endl;
     like << " " << "pen_prop_aco" << " " << pen_prop_aco << endl;
     like << " " << "pen_cb_sel" << " " << pen_cb_sel << endl;
     like << " " << "pen_rdev" << " " << pen_rdev << endl;
     like << " " << "pen_fdev" << " " << pen_fdev << endl;
     //like << " " << "total_survey_index" << " " << Lindex_md+sum(Lindex_cm)+Lindex_ny+Lindex_nj+sum(Lindex_ct)+Lindex_des+Lindex_de30 << endl;
     //like << " " << "total_survey_agecomp" << " " << Lindexagecomp_md+sum(Lindexagecomp_cm)+Lindexagecomp_ny+Lindexagecomp_nj+sum(Lindexagecomp_ct)+ Lindexagecomp_des+Lindexagecomp_de30 << endl;
     //like << " " << "total_age1_index" << " " << Lage1index_md+sum(Lage1index_ny) << endl;
     //like << " " << "total_yoy_index" << " " << sum(Lyoyindex_bay)+sum(Lyoyindex_coast) << endl;
     //like << " " << "total_neg_ll" << " " << neg_LL << endl;
   }
   ofstream occprob("occ.txt");//occupancy probability
   {
     occprob << "param stock timestep age prob logsd" << endl;
       //for(s=1;s<=stock;s++)
       //{
         for(t=1;t<=tstep;t++)
         {
           for(a=fage+1;a<=lage;a++)
           {
              occprob << " " << "est" <<  " " << "1" << " " << t << " " << a  << " " << log_prop_bay(t,a) << " " << log_sd_bay(t,a) << endl;
              occprob << " " << "obs" << " " << "1" << " " << t << " " << a << " " <<  prop_bay(t,a) << " " << log_sd_bay(t,a) << endl;
              occprob << " " << "est" << " " << "2" << " " << t << " " << a  << " " << log_prop_coast(t,a) << " " << log_sd_coast(t,a) << endl;
              occprob << " " << "obs" << " " << "2" << " " << t << " " << a << " " <<  prop_coast(t,a) << " " << log_sd_coast(t,a) << endl;
           }//close age+1 loop
         }//close tstep 
       //}//close stock loop
    }//close ofstream
   ofstream ssb("ssb.txt");//spawning stock biomass
   {
     ssb << "stock region year ssb" << endl;
       for(s=1;s<=stock;s++)
       {
         for(r=1;r<=region;r++)
         {
           for(y=fmyear;y<=lmyear;y++)
           {
              ssb << " " << s <<  " " << r << " " << y  << " " << SSB(s,r,y) << endl;
           }//close year loop
         }//close region 
       }//close stock loop
    }//close ofstream
   ofstream bio("biomass.txt");//spawning stock biomass
   {
     bio << "stock region year bio" << endl;
       for(s=1;s<=stock;s++)
       {
         for(r=1;r<=region;r++)
         {
           for(y=fmyear;y<=lmyear;y++)
           {
             //for(ts=1;ts<=tstep;ts++)
             //{
              bio << " " << s <<  " " << r << " " << y  << " " << B(s,r,y) << endl;
             //}
           }//close year loop
         }//close region 
       }//close stock loop
    }//close ofstream
}

void model_parameters::set_runtime(void)
{
  dvector temp1("{20000, 20000, 20000 //change the maximum number of iterations for each phase}");
  maximum_function_evaluations.allocate(temp1.indexmin(),temp1.indexmax());
  maximum_function_evaluations=temp1;
}

void model_parameters::preliminary_calculations(void){
#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::final_calcs(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
  gradient_structure::set_MAX_NVAR_OFFSET(1000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(100000);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(100000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(1000000);
  arrmblsize=900000;
    gradient_structure::set_NO_DERIVATIVES();
#ifdef DEBUG
  #ifndef __SUNPRO_C
std::feclearexcept(FE_ALL_EXCEPT);
  #endif
#endif
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
#ifdef DEBUG
  #ifndef __SUNPRO_C
bool failedtest = false;
if (std::fetestexcept(FE_DIVBYZERO))
  { cerr << "Error: Detected division by zero." << endl; failedtest = true; }
if (std::fetestexcept(FE_INVALID))
  { cerr << "Error: Detected invalid argument." << endl; failedtest = true; }
if (std::fetestexcept(FE_OVERFLOW))
  { cerr << "Error: Detected overflow." << endl; failedtest = true; }
if (std::fetestexcept(FE_UNDERFLOW))
  { cerr << "Error: Detected underflow." << endl; }
if (failedtest) { std::abort(); } 
  #endif
#endif
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
