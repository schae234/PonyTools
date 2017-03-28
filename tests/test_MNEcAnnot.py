import ponytools as pc
annot = pc.MNEc2MAnnot()

def test_init():
    assert all([x.startswith('AX') for x in annot._annot.ProbeSetID])
    assert all([x.startswith('Affx') for x in annot._annot.AffySNPID])
    assert all([x.startswith('MNEc') for x in annot._annot.MNEcID])
