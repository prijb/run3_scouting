cmsRun producer_Run3.py inputs=/ceph/cms/store/user/mmasciov/nfs-7/ScoutingPFRun3_Run2022D_forTest.root era=2022D data=True nevents=1000 output=output_data.root >& log_data.txt&
cmsRun producer_Run3.py inputs=/ceph/cms/store/user/mmasciov/ScoutingPFRun3_signalMCScenarioA_default_forTest.root era=2022 data=False nevents=1000 output=output_mc.root >& log_mc.txt&

cmsRun producer_Run3.py inputs=/ceph/cms/store/user/legianni/test2023-scouting/Run2023B-379d9f7f-67f2-46c9-987b-af690532b4f9.root era=2023B data=True nevents=1000 output=output_data_run3B.root >& log_data.txt&
cmsRun producer_Run3.py inputs=/ceph/cms/store/user/legianni/test2023-scouting/Run2023C-367094-0e43dc5a-3fd6-439e-9f9b-73bec29a56d9.root era=2023C-triggerV10 data=True nevents=1000 output=output_data_run3CB.root >& log_data.txt&
cmsRun producer_Run3.py inputs=/ceph/cms/store/user/legianni/test2023-scouting/Run2023C-367838-00a80575-133a-4d2c-a2ac-a8ab44f71074.root era=2023C data=True nevents=1000 output=output_data_run3C.root >& log_data.txt&


