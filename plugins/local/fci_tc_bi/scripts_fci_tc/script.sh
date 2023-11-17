source /home_lct/eginer/qp2/quantum_package.rc
input=$1
 basis=$2
 atom=$3
 mul=$4
 EXPORT_OMP_NUM_THREADS=16
 dir=${input}_${basis}
 mkdir ${dir}
 cp ${input}.xyz ${dir}/
 cd $dir
 EZFIO=${input}_${basis}_bi_ortho
 qp create_ezfio -b "${atom}:cc-pcv${basis}|H:cc-pv${basis}" ${input}.xyz -m $mul -o $EZFIO
 qp run scf
 # Getting THE GOOD VALUE OF MU
 qp run print_mu_av_tc  | tee ${EZFIO_FILE}.mu_av.out
 mu=`grep "average_mu_rs_c_lda  =" ${EZFIO_FILE}.mu_av.out | cut -d "=" -f 2`
 qp set ao_two_e_erf_ints mu_erf $mu
 # Carrying the BI-ORTHO TC-SCF 
 qp run tc_scf | tee ${EZFIO_FILE}.tc_scf.out
 #   Three body terms without normal order
 ### THREE E TERMS FOR FCI
 qp set tc_keywords three_body_h_tc True  
 qp set tc_keywords double_normal_ord False
 qp set perturbation pt2_max 0.003
 qp run fci_tc_bi_ortho | tee ${EZFIO_FILE}.tc_fci.out
 #   Three body terms with normal order
 qp set tc_keywords double_normal_ord True
 qp run fci_tc_bi_ortho | tee ${EZFIO_FILE}.tc_fci_normal_order.out

cd ../

