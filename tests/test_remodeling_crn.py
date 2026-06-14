from afib_microgravity.remodeling import baseline_crn_params, microgravity_crn_params


def test_af_type_remodelling_scales_intended_currents_down():
    base = baseline_crn_params()
    mg = microgravity_crn_params(severity=1.0)
    assert mg.g_CaL < base.g_CaL          # ICaL down
    assert mg.g_to  < base.g_to           # Ito down
    assert mg.g_Kur_scale < base.g_Kur_scale  # IKur down (scale multiplier)
    assert mg.g_Na  == base.g_Na          # INa untouched
    assert mg.g_Kr  == base.g_Kr          # IKr untouched


def test_severity_zero_recovers_baseline():
    base = baseline_crn_params()
    mg   = microgravity_crn_params(severity=0.0)
    assert mg.g_CaL       == base.g_CaL
    assert mg.g_to        == base.g_to
    assert mg.g_Kur_scale == base.g_Kur_scale
