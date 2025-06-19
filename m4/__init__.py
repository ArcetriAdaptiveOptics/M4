
from opticalib import PhaseCam
from m4.configuration import ott
interf = PhaseCam('6110')
ott = ott.create_ott()


try:
    from opticalib import AdOpticaDm
    dm = AdOpticaDm()
except Exception as e:
    pass
    # print(e)
    # print("simulating M4 DM")
    # from m4.simulator import fake_deformable_mirror
    # dm = fake_deformable_mirror.FakeM4DM()

