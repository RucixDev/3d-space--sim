import os
import random
import unittest

from tools.stellar_content_pipeline.schema import load_universe_state_json, universe_state_to_dict, UniverseState
from tools.stellar_content_pipeline.risk import RiskModel, RoutePlanner
from tools.stellar_content_pipeline.galnet import GalNetSynth
from tools.stellar_content_pipeline.missions import MissionSynth
from tools.stellar_content_pipeline.simulator import TimelineSimulator


class PipelineTests(unittest.TestCase):
    def setUp(self) -> None:
        root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
        self.example = os.path.join(root, "tools", "examples", "universe_state.json")
        self.state = load_universe_state_json(self.example)

    def test_deterministic_bulletins(self) -> None:
        g1 = GalNetSynth(seed=42)
        g2 = GalNetSynth(seed=42)
        b1 = g1.generate_bulletins(self.state, max_items=10)
        b2 = g2.generate_bulletins(self.state, max_items=10)
        self.assertEqual([x.id for x in b1], [x.id for x in b2])

    def test_route_planning(self) -> None:
        rp = RoutePlanner()
        # Aquila -> Juno should be reachable in example
        route = rp.plan(self.state, 10, 19)
        self.assertIsNotNone(route)
        assert route is not None
        self.assertEqual(route[0], 10)
        self.assertEqual(route[-1], 19)

    def test_mission_generation(self) -> None:
        rm = RiskModel()
        g = GalNetSynth(seed=7, risk_model=rm)
        bulletins = g.generate_bulletins(self.state, max_items=8)

        ms = MissionSynth(seed=7, risk_model=rm)
        missions = ms.generate(self.state, bulletins, max_items=12)
        self.assertGreaterEqual(len(missions), 5)
        # Ensure essential fields look sane
        for m in missions:
            self.assertTrue(m.id)
            self.assertGreater(m.reward_credits, 0)
            self.assertGreaterEqual(m.deadline_day, self.state.day)

    def test_timeline_simulation(self) -> None:
        sim = TimelineSimulator(seed=123)
        timeline = sim.simulate(self.state, days=3, step_days=1.0, bulletins_per_day=6, missions_per_day=6)
        self.assertEqual(len(timeline), 3)
        # Days should strictly increase
        days = [st.day for (st, _, _) in timeline]
        self.assertTrue(days[0] < days[1] < days[2])

    def test_universe_round_trip(self) -> None:
        d = universe_state_to_dict(self.state)
        st2 = UniverseState.from_json(d)
        self.assertAlmostEqual(st2.day, self.state.day, places=6)
        self.assertEqual([s.id for s in st2.systems], [s.id for s in self.state.systems])
        self.assertEqual([f.id for f in st2.factions], [f.id for f in self.state.factions])

    def test_advance_state(self) -> None:
        sim = TimelineSimulator(seed=123)
        rng = random.Random(999)
        st2 = sim.advance_state(self.state, days=2.5, step=1.0, rng=rng)
        self.assertAlmostEqual(st2.day, self.state.day + 2.5, places=6)


if __name__ == "__main__":
    unittest.main()
