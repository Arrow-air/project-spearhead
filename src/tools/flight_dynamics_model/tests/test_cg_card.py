import pytest

from spearhead.cg_card import cg_card, cg_card_text


def test_cg_card_matches_the_window_analysis_at_the_target_cg():
    card = cg_card("pt15", -.04, .065)
    assert card["static_margin_22"] == pytest.approx(.064, abs=.004)
    # 22 kg with the MAD lift units and HOBBYWING pusher (note 0003 at 21 kg: -18.2°).
    assert card["trims"][22.]["elevator_deg"] == pytest.approx(-18.4, abs=.2)
    assert card["trims"][22.]["beyond_measured"] and card["trims"][24.]["beyond_measured"]
    assert .45 < card["hover_max_rotor_fraction"] < .5
    assert card["pusher"]["option"] == "hobbywing"
    assert .7 < card["trims"][22.]["pusher_command"] < .8 and not card["trims"][22.]["above_current_limit"]
    assert card["rotor_centroid_x_m"] == pytest.approx(-.0375)
    text = cg_card_text(card)
    assert "SM at 22 m/s: +6.4 %" in text and "beyond the measured ±15°" in text


def test_moving_cg_aft_lowers_margin_and_elevator_demand():
    aft, target = cg_card("pt15", -.06, .065), cg_card("pt15", -.04, .065)
    assert aft["static_margin_22"] < target["static_margin_22"]
    assert aft["trims"][22.]["elevator_deg"] > target["trims"][22.]["elevator_deg"]
