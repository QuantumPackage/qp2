QP_ROOT=/home/eginer/new_qp2/qp2
source ${QP_ROOT}/quantum_package.rc
  echo Ne > Ne.xyz
  echo $QP_ROOT
  qp create_ezfio -b cc-pcvdz Ne.xyz 
  qp run scf 
  qp set tc_keywords bi_ortho True 
  qp set ao_two_e_erf_ints mu_erf 0.87 
  qp set tc_keywords j1b_pen [1.5]
  qp set tc_keywords j1b_type 3 
  qp run tc_scf | tee Ne.ezfio.tc_scf.out 
  grep "TC energy =" Ne.ezfio.tc_scf.out | tail -1 
  eref=-128.552134
