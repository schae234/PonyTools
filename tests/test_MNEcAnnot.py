import ponytools as pc
annot = pc.MNEc2MAnnot()

def test_init():
    assert all([x.startswith('AX') for x in annot._annot.reset_index().ProbeSetID])
    assert all([x.startswith('Affx') for x in annot._annot.reset_index().AffySNPID])
    assert all([x.startswith('MNEc') for x in annot._annot.reset_index().MNEcID])


def test_in2M():
    assert annot.in2M('chr23',22999655)
    assert annot.in2M('chr23',22999656) == False

    assert annot.in2M('chr18',66493737)
    assert annot.in2M('chr18',66493738) == False


def test_in670():
    assert annot.in670('chr23',22999655) == False

    assert annot.in670('chr18',66493737)
