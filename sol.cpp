#include <iostream>

#include "dbg.h"
#include "json.hpp"

using namespace std;

struct Point {
	double x, y;
};

struct SubstationType {
	double cs, rs, ps;
};

struct CableType {
	double cqf, cql, pq, rq;
};

struct Scenario {
	double pw, piw;
};

struct Answer {
	struct Station {
		int id, st=-1, lct, other=-1, sst;
	};
	vector<Station> stations;
	vector<int> turbines;
};

double ctf, ctl;
double c0, cP, Cmax;
double Pmax;
Point v0;
vector<Point> substations, turbines;
vector<SubstationType> sts;
vector<CableType> ls_cts, ss_cts;
vector<Scenario> scenarios;
int S, T;

void read() {
	nlohmann::json js;
	cin >> js;

	{const auto &gp = js["general_parameters"];
		c0 = gp["curtailing_cost"];
		cP = gp["curtailing_penalty"];
		Cmax = gp["maximum_curtailing"];
		ctf = gp["fixed_cost_cable"];
		ctl = gp["variable_cost_cable"];
		v0.x = gp["main_land_station"]["x"];
		v0.y = gp["main_land_station"]["y"];
		Pmax = gp["maximum_power"];
	}

	substations.resize(S = js["substation_locations"].size());
	for(const auto &sl : js["substation_locations"]) {
		Point &s = substations[int(sl["id"])-1];
		s.x = sl["x"];
		s.y = sl["y"];
	}

	turbines.resize(T = js["wind_turbines"].size());
	for(const auto &wt0 : js["wind_turbines"]) {
		Point &t = turbines[int(wt0["id"])-1];
		t.x = wt0["x"];
		t.y = wt0["y"];
	}

	sts.resize(js["substation_types"].size());
	for(const auto &st0 : js["substation_types"]) {
		SubstationType &st = sts[int(st0["id"])-1];
		st.cs = st0["cost"];
		st.ps = st0["probability_of_failure"];
		st.rs = st0["rating"];
	}

	ls_cts.resize(js["land_substation_cable_types"].size());
	for(const auto &ct0 : js["land_substation_cable_types"]) {
		CableType &ct = ls_cts[int(ct0["id"])-1];
		ct.cqf = ct0["fixed_cost"];
		ct.cql = ct0["variable_cost"];
		ct.pq = ct0["probability_of_failure"];
		ct.rq = ct0["rating"];
	}

	ss_cts.resize(js["substation_substation_cable_types"].size());
	for(const auto &ct0 : js["substation_substation_cable_types"]) {
		CableType &ct = ss_cts[int(ct0["id"])-1];
		ct.cqf = ct0["fixed_cost"];
		ct.cql = ct0["variable_cost"];
		ct.rq = ct0["rating"];
	}

	scenarios.resize(js["wind_scenarios"].size());
	for(const auto &ws0 : js["wind_scenarios"]) {
		Scenario &s = scenarios[int(ws0["id"])-1];
		s.pw = ws0["probability"];
		s.piw = ws0["power_generation"];
	}
}

void write(const Answer &ans) {
	nlohmann::json js;
	js["substation_substation_cables"] = nlohmann::json::array();

	for(const auto &s : ans.stations) {
		if(s.st < 0) continue;
		nlohmann::json j;
		j["id"] = s.id+1;
		j["land_cable_type"] = s.lct+1;
		j["substation_type"] = s.st+1;
		js["substations"].push_back(j);
		if(s.other > s.id) {
			nlohmann::json j;
			j["substation_id"] = s.id+1;
			j["other_substation_id"] = s.other+1;
			j["cable_type"] = s.sst+1;
			js["substation_substation_cables"].push_back(j);
		}
	}

	for(int i = 0; i < T; ++i) {
		nlohmann::json j;
		j["id"] = i+1;
		j["substation_id"] = ans.turbines[i]+1;
		js["turbines"].push_back(j);
	}

	cout << js << endl;
}

double cC(double C) {
	double ans = c0 * C;
	if(C > Cmax) ans += cP * (C - Cmax);
	return ans; 
}

double dist(const Point &a, const Point &b) {
	const double dx = b.x - a.x;
	const double dy = b.y - a.y;
	return sqrt(dx*dx+dy*dy);
}

double ceq(const Point &a, const Point &b, const CableType &q) {
	return q.cqf + q.cql * dist(a, b);
}

double ce(int i, int j) {
	return ctf + ctl * dist(substations[i], turbines[j]);
}

double pf(const Answer &ans) {
	double p = 0.;
	for(const auto &s : ans.stations) if(s.st >= 0) {
		p += sts[s.st].ps + ls_cts[s.lct].pq;
	}
	return p;
}

double Cn(const Answer &ans, double pi) {
	double C = 0.;
	vector<double> pv(S, 0.);
	for(int s : ans.turbines) pv[s] += pi;
	for(const auto &s : ans.stations) if(s.st >= 0) {
		const double m = min(sts[s.st].rs, ls_cts[s.lct].rq);
		if(pv[s.id] > m) C += pv[s.id] - m;
	}
	return C;
}

double Cf(const Answer &ans, double pi) {
	double C = 0.;
	vector<double> pv(S, 0.);
	for(int s : ans.turbines) pv[s] += pi;
	for(const auto &s : ans.stations) if(s.st >= 0) {
		const double pf = sts[s.st].ps + ls_cts[s.lct].pq;
		if(s.other < 0) {
			C += pf * cC(pv[s.id]);
		} else {
			double cf = 0.;
			const double r = ss_cts[s.sst].rq;
			if(pv[s.id] > r) cf += pv[s.id] - r;
			const double x = pv[s.other] + min(r, pv[s.id]) - min(sts[ans.stations[s.other].st].rs, ls_cts[ans.stations[s.other].lct].rq);
			if(x > 0) cf += x;
			C += pf * cC(cf);
		}
	}
	return C;
}

bool show = false;
double c(const Answer &ans) {
	double c = 0.;

	for(const auto &s : ans.stations) if(s.st >= 0) {
		c += sts[s.st].cs + ceq(v0, substations[s.id], ls_cts[s.lct]);
		if(s.other > s.id) c += ceq(substations[s.id], substations[s.other], ss_cts[s.sst]);
	}
	if(show) dbg(c);
	for(int i = 0; i < T; ++i) {
		c += ce(ans.turbines[i], i);
	}
	if(show) dbg(c);

	const double pfa = pf(ans);
	for(const Scenario &s : scenarios) {
		c += s.pw * (Cf(ans, s.piw) + (1-pfa) * cC(Cn(ans, s.piw)));
	}

	return c;
}

int main() {
	ios::sync_with_stdio(false);
	cin.tie(nullptr);
	read();

	dbg(S, T, sts.size(), ls_cts.size(), ss_cts.size(), scenarios.size(), ctf, ctl, c0, cP);

	Point mean = {0., 0.};
	for(const Point &p : turbines) {
		mean.x += p.x;
		mean.y += p.y;
	}
	mean.x /= T;
	mean.y /= T;

	Answer best;
	double score = 1e9;

	for(int st = 0; st < sts.size(); ++st) for(int lct = 0; lct < ls_cts.size(); ++lct) {
		if(sts[st].rs != ls_cts[lct].rq) continue;
		const int SMax = clamp(S, 1, min(S, (int) ceil(Pmax * T / sts[st].rs)));
		for(int S2 = 1; S2 <= SMax; ++S2) {
			Answer ans;
			ans.stations.resize(S);
			ans.turbines.resize(T);
			for(int i = 0; i < S; ++i) ans.stations[i].id = i;
			vector<pair<double, int>> bs;
			for(int i = 0; i < S; ++i) {
				const double d = S2 * ceq(v0, substations[i], ls_cts[lct]) + T * (ctf + ctl * dist(substations[i], mean));
				bs.emplace_back(d, i);
			}
			ranges::sort(bs);

			for(int i = 0; i < S2; ++i) {
				ans.stations[bs[i].second].st = st;
				ans.stations[bs[i].second].lct = lct;
			}

			vector<tuple<double, int, int>> match;
			for(int i = 0; i < S2; ++i) for(int j = 0; j < T; ++j) {
				const double d = ce(bs[i].second, j);
				match.emplace_back(d, bs[i].second, j);
			}
			ranges::sort(match);
			ans.turbines.assign(T, -1);
			vector<int> count(S, 0);
			int left = T%S2, mt = T/S2;
			for(auto &[d, i, j] : match) {
				if(ans.turbines[j] >= 0) continue;
				if(count[i] > mt) continue;
				if(count[i] == mt) {
					if(left) --left;
					else continue;
				}
				++ count[i];
				ans.turbines[j] = i;
			}

			double sc = c(ans);
			if(sc < score) {
				if(sc < 9000) dbg(st, lct, sts[st].rs, sc, S2, SMax, mt*Pmax, (mt+1)*Pmax);
				best = ans;
				score = sc;
			}

			vector<tuple<double, int, int>> ssd;
			for(int i = 0; i < S2; ++i) for(int j = i+1; j < S2; ++j) {
				ssd.emplace_back(dist(substations[bs[i].second], substations[bs[j].second]),
									bs[i].second,
									bs[j].second);
			}
			ranges::sort(ssd);
			for(auto &[d, i, j] : ssd) {
				if(ans.stations[i].other >= 0) continue;
				if(ans.stations[j].other >= 0) continue;
				ans.stations[i].other = j;
				ans.stations[j].other = i;
			}

			for(int sst = 0; sst < ss_cts.size(); ++sst) {
				for(auto &s : ans.stations) if(s.other >= 0) s.sst = sst;
				sc = c(ans);
				if(sc < score) {
					if(sc < 9000) dbg(st, lct, sts[st].rs, sc, S2, SMax, mt*Pmax, (mt+1)*Pmax, "HERE", ss_cts[sst].rq);
					best = ans;
					score = sc;
				}
			}
		}
	}	

	show = true;
	dbg(c(best));
	write(best);

	return 0;
}