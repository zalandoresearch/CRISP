T = 150
S = 1000
P0 = 0.0001
P1 = 0.025
TEST_CAPACITY = 0.01

R0 = 2.5 # 5
POLICY_START = 30 # 60 
IT = 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
SEED = 42 # 4711 0815 666 12345 13156 80808 1000 0 1024 # 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
PEIT = 0.5
PSRT = 0.6 0.7 0.8
QDAYS = 2 7 14 21
PDAYS = 7

ARGS:= --R0_{R0}_--interactive_off_--policy-start_{PS}_--p0_$(P0)_--p1_$(P1)_--test-capacity_$(TEST_CAPACITY)_--T_$(T)_--S_$(S)

BUILD_TARGETS_NO:= _--policy_no_--seed_{SD}_--iteration_{IT}_${ARGS}
BUILD_TARGETS_NO:=$(foreach IT,$(IT),$(subst {IT},$(IT),$(BUILD_TARGETS_NO)))
BUILD_TARGETS_NO:=$(foreach R0,$(R0),$(subst {R0},$(R0),$(BUILD_TARGETS_NO)))
BUILD_TARGETS_NO:=$(foreach PS,$(POLICY_START),$(subst {PS},$(PS),$(BUILD_TARGETS_NO)))
BUILD_TARGETS_NO:=$(foreach SD,$(SEED), $(subst {SD},$(SD),$(BUILD_TARGETS_NO)))

BUILD_TARGETS_LOCKDOWN:= _--policy_lockdown_--seed_{SD}_--iteration_{IT}_${ARGS}
BUILD_TARGETS_LOCKDOWN:=$(foreach IT,$(IT),$(subst {IT},$(IT),$(BUILD_TARGETS_LOCKDOWN)))
BUILD_TARGETS_LOCKDOWN:=$(foreach R0,$(R0),$(subst {R0},$(R0),$(BUILD_TARGETS_LOCKDOWN)))
BUILD_TARGETS_LOCKDOWN:=$(foreach PS,$(POLICY_START),$(subst {PS},$(PS),$(BUILD_TARGETS_LOCKDOWN)))
BUILD_TARGETS_LOCKDOWN:=$(foreach SD,$(SEED), $(subst {SD},$(SD),$(BUILD_TARGETS_LOCKDOWN)))

BUILD_TARGETS_SYMPTOM:= _--policy_symptom_--qDays_{QD}_--seed_{SD}_--iteration_{IT}_${ARGS}
BUILD_TARGETS_SYMPTOM:=$(foreach IT,$(IT),$(subst {IT},$(IT),$(BUILD_TARGETS_SYMPTOM)))
BUILD_TARGETS_SYMPTOM:=$(foreach R0,$(R0),$(subst {R0},$(R0),$(BUILD_TARGETS_SYMPTOM)))
BUILD_TARGETS_SYMPTOM:=$(foreach PS,$(POLICY_START),$(subst {PS},$(PS),$(BUILD_TARGETS_SYMPTOM)))
BUILD_TARGETS_SYMPTOM:=$(foreach QD,$(QDAYS),$(subst {QD},$(QD),$(BUILD_TARGETS_SYMPTOM)))
BUILD_TARGETS_SYMPTOM:=$(foreach SD,$(SEED), $(subst {SD},$(SD),$(BUILD_TARGETS_SYMPTOM)))

BUILD_TARGETS_CONTACT:= _--policy_contact_--qDays_{QD}_--pDays_{PD}_--seed_{SD}_--iteration_{IT}_${ARGS}
BUILD_TARGETS_CONTACT:=$(foreach IT,$(IT),$(subst {IT},$(IT),$(BUILD_TARGETS_CONTACT)))
BUILD_TARGETS_CONTACT:=$(foreach R0,$(R0),$(subst {R0},$(R0),$(BUILD_TARGETS_CONTACT)))
BUILD_TARGETS_CONTACT:=$(foreach PS,$(POLICY_START),$(subst {PS},$(PS),$(BUILD_TARGETS_CONTACT)))
BUILD_TARGETS_CONTACT:=$(foreach QD,$(QDAYS),$(subst {QD},$(QD),$(BUILD_TARGETS_CONTACT)))
BUILD_TARGETS_CONTACT:=$(foreach PD,$(PDAYS),$(subst {PD},$(PD),$(BUILD_TARGETS_CONTACT)))
BUILD_TARGETS_CONTACT:=$(foreach SD,$(SEED), $(subst {SD},$(SD),$(BUILD_TARGETS_CONTACT)))

BUILD_TARGETS_SCORE:= _--policy_lbp_--pEI-threshold_{PEIT}_--pSR-threshold_{PSRT}_--seed_{SD}_--iteration_{IT}_${ARGS}
BUILD_TARGETS_SCORE:=$(foreach IT,$(IT),$(subst {IT},$(IT),$(BUILD_TARGETS_SCORE)))
BUILD_TARGETS_SCORE:=$(foreach R0,$(R0),$(subst {R0},$(R0),$(BUILD_TARGETS_SCORE)))
BUILD_TARGETS_SCORE:=$(foreach PS,$(POLICY_START),$(subst {PS},$(PS),$(BUILD_TARGETS_SCORE)))
BUILD_TARGETS_SCORE:=$(foreach PEIT,$(PEIT),$(subst {PEIT},$(PEIT),$(BUILD_TARGETS_SCORE)))
BUILD_TARGETS_SCORE:=$(foreach PSRT,$(PSRT),$(subst {PSRT},$(PSRT),$(BUILD_TARGETS_SCORE)))
BUILD_TARGETS_SCORE:=$(foreach SD,$(SEED), $(subst {SD},$(SD),$(BUILD_TARGETS_SCORE)))

# BUILD_TARGETS:=$(BUILD_TARGETS_NO) $(BUILD_TARGETS_LOCKDOWN) $(BUILD_TARGETS_SYMPTOM) $(BUILD_TARGETS_CONTACT) 
BUILD_TARGETS:=$(BUILD_TARGETS_NO) $(BUILD_TARGETS_LOCKDOWN) $(BUILD_TARGETS_SYMPTOM) $(BUILD_TARGETS_CONTACT) $(BUILD_TARGETS_SCORE)

# uncomment this one when not running on SLURM
SRUN := srun -p cpu -c 1 --mem=8GB
# SRUN := 

CMD := python -u ../../exp_test_quarantine.py

RUN_EXP = run_exp() { \
	mkdir results/$$1; \
	ln -sf results/$$1 $$2.part; \
	cd results/$$1; \
	${SRUN} ${CMD} $$(echo "$$2" | sed -e 's/_/ /g') >sout 2>serr  && \
	cd ../.. && \
	mv $$2.part $$2; \
}

ifneq (,$(findstring n,$(MAKEFLAGS)))
RUN_EXP=: RUN_EXP
endif

tasks: $(BUILD_TARGETS)
.PHONY: clean tasks dist-clean


$(BUILD_TARGETS): TMP := $(shell mktemp -u XXXXXXXXXX)

$(BUILD_TARGETS):
	$(RUN_EXP); run_exp $(TMP) $@


results:
	mkdir results

clean:
	for i in *.part; do  rm -r $$(readlink $$i);  rm $$i; done

dist-clean:
	rm -rf results/* _--*
