from afib_microgravity.restitution import apd_restitution


def test_apd_restitution_is_monotonic_and_shortens_with_short_di():
    di, apd = apd_restitution(diastolic_intervals=(300., 200., 120., 60.), dt=0.05)
    assert all(a > 0 for a in apd)
    assert apd[0] > apd[-1]                                   # shorter DI -> shorter APD
    assert all(apd[i] >= apd[i+1] for i in range(len(apd)-1)) # monotonic non-increasing
