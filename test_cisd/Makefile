FC      = gfortran
FFLAGS  = -O2 -Wall
LDFLAGS = -L$(HOME)/miniconda3/envs/libra-devel/lib -lopenblas -Wl,-rpath,$(HOME)/miniconda3/envs/libra-devel/lib
MOD_DIR = .
OBJS_B  = cis_overlap_cisd.o cis_overlap_cisd_l4m.o cis_overlap_cisd_l4m_fast.o benchmark_l4m.o

benchmark_l4m: $(OBJS_B)
	$(FC) $(FFLAGS) -o $@ $(OBJS_B) $(LDFLAGS)

cis_overlap_cisd.o: ../src/util_chem/cis/cis_overlap_cisd.f90
	$(FC) $(FFLAGS) -J$(MOD_DIR) -c $< -o $@

cis_overlap_cisd_l4m.o: ../src/util_chem/cis/cis_overlap_cisd_l4m.f90 cis_overlap_cisd.o
	$(FC) $(FFLAGS) -I$(MOD_DIR) -J$(MOD_DIR) -c $< -o $@

cis_overlap_cisd_l4m_fast.o: ../src/util_chem/cis/cis_overlap_cisd_l4m_fast.f90 cis_overlap_cisd.o
	$(FC) $(FFLAGS) -I$(MOD_DIR) -J$(MOD_DIR) -c $< -o $@

benchmark_l4m.o: benchmark_l4m.f90
	$(FC) $(FFLAGS) -I$(MOD_DIR) -c $< -o $@

clean:
	rm -f benchmark_l4m *.o *.mod

.PHONY: all clean

check_l4m: check_l4m.o cis_overlap_cisd.o cis_overlap_cisd_l4m.o cis_overlap_cisd_l4m_fast.o
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)
check_l4m.o: check_l4m.f90
	$(FC) $(FFLAGS) -I$(MOD_DIR) -c $< -o $@

check_l4m_extended: check_l4m_extended.o cis_overlap_cisd.o cis_overlap_cisd_l4m.o cis_overlap_cisd_l4m_fast.o
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)
check_l4m_extended.o: check_l4m_extended.f90
	$(FC) $(FFLAGS) -I$(MOD_DIR) -c $< -o $@

check_real: check_real.o cis_overlap_cisd.o cis_overlap_cisd_l4m.o cis_overlap_cisd_l4m_fast.o
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)
check_real.o: check_real.f90
	$(FC) $(FFLAGS) -I$(MOD_DIR) -c $< -o $@

check_multi: check_multi.o cis_overlap_cisd.o cis_overlap_cisd_l4m.o cis_overlap_cisd_l4m_fast.o
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)
check_multi.o: check_multi.f90
	$(FC) $(FFLAGS) -I$(MOD_DIR) -c $< -o $@

check_trunc: check_trunc.o cis_overlap_cisd.o cis_overlap_cisd_l4m.o cis_overlap_cisd_l4m_fast.o
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)
check_trunc.o: check_trunc.f90
	$(FC) $(FFLAGS) -I$(MOD_DIR) -c $< -o $@

check_real_tri: check_real_tri.o cis_overlap_cisd.o cis_overlap_cisd_l4m.o cis_overlap_cisd_l4m_fast.o
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)
check_real_tri.o: check_real_tri.f90
	$(FC) $(FFLAGS) -I$(MOD_DIR) -c $< -o $@
