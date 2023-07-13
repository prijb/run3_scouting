cmsRun producer_Run3.py inputs=/ceph/cms/store/user/mmasciov/nfs-7/ScoutingPFRun3_Run2022D_forTest.root era=2022D data=True nevents=1000 output=output_data.root >& log_data.txt&
cmsRun producer_Run3.py inputs=/ceph/cms/store/user/mmasciov/ScoutingPFRun3_signalMCScenarioA_default_forTest.root era=2022 data=False nevents=1000 output=output_mc.root >& log_mc.txt&
