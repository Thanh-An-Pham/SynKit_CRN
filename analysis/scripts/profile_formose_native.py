"""Instrument a private source copy and profile native formose generation.

Production sources and the paper's benchmark reports are left untouched.
Timers are inclusive and profiling is restricted to one worker thread.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import difflib
import hashlib
import json
import os
from pathlib import Path
import shutil
import statistics
import sys
import time


def instrument(source):
    def replace(old,new):
        nonlocal source
        assert source.count(old)==1,(old[:100],source.count(old))
        source=source.replace(old,new)
    replace('class Expander {', '''struct ProfileEntry { double seconds = 0; std::int64_t calls = 0; };
std::map<std::string, ProfileEntry> profile_entries;
struct ProfileScope {
  const char *name;
  Clock::time_point start = Clock::now();
  explicit ProfileScope(const char *name) : name(name) {}
  ~ProfileScope() {
    auto &entry = profile_entries[name];
    entry.seconds += std::chrono::duration<double>(Clock::now() - start).count();
    ++entry.calls;
  }
};

class Expander {''')
    # Inclusive totals, split further using the scopes below.
    replace('    auto cached = standard_.find(text);','    ProfileScope scope("standardize");\n    auto cached = standard_.find(text);')
    replace('    Species record;','    ProfileScope scope("new_species");\n    Species record;')
    replace('    for (std::size_t r = 0; r < rules_.size(); ++r) {','    ProfileScope binding_scope("new_species_bindings");\n    for (std::size_t r = 0; r < rules_.size(); ++r) {')
    replace('    auto mol = RDKit::v2::SmilesParse::MolFromSmiles(key);','    auto mol = [&]() { ProfileScope s("new_species_parse"); return RDKit::v2::SmilesParse::MolFromSmiles(key); }();')
    replace('    RDKit::MolOps::addHs(*explicit_mol);','    { ProfileScope s("new_species_add_hs"); RDKit::MolOps::addHs(*explicit_mol); }')
    replace('    std::vector<int> ranked;','    ProfileScope scope("make_tasks");\n    std::vector<int> ranked;')
    replace('    TaskResult result;','    ProfileScope scope("apply");\n    TaskResult result;')
    replace('      auto outcomes = rule.active().runReactants(reactants, o_.max_products);','      auto outcomes = [&]() { ProfileScope s("run_reactants"); return rule.active().runReactants(reactants, o_.max_products); }();')
    replace('    std::vector<std::string> parts;','    ProfileScope scope("canonical_outcome");\n    std::vector<std::string> parts;')
    replace('      RDKit::MolOps::sanitizeMol(molecule);','      { ProfileScope s("outcome_sanitize"); RDKit::MolOps::sanitizeMol(molecule); }')
    replace('      if (!o_.explicit_h) RDKit::MolOps::removeHs(molecule, RDKit::MolOps::RemoveHsParameters(), true);','      if (!o_.explicit_h) { ProfileScope s("outcome_remove_hs"); RDKit::MolOps::removeHs(molecule, RDKit::MolOps::RemoveHsParameters(), true); }')
    replace('      parts.push_back(RDKit::MolToSmiles(molecule));','      { ProfileScope s("outcome_smiles"); parts.push_back(RDKit::MolToSmiles(molecule)); }')
    replace('    std::vector<TaskResult> results(tasks.size());','    ProfileScope scope("run_tasks");\n    std::vector<TaskResult> results(tasks.size());')
    replace('    std::set<std::string> next;','    ProfileScope scope("integrate");\n    std::set<std::string> next;')
    replace('    json::object out;\n    out["protocol"] = kProtocol;','    json::object out;\n    out["protocol"] = kProtocol;') if False else None
    replace('    out["rules"] = std::move(rules);\n    return out;','''    out["rules"] = std::move(rules);
    json::object profile;
    for (const auto &[name, entry] : profile_entries) {
      json::object item;
      item["seconds"] = entry.seconds;
      item["calls"] = entry.calls;
      profile[name] = std::move(item);
    }
    out["profile"] = std::move(profile);
    return out;''')
    replace('    const json::object request = json::parse(text).as_object();','    profile_entries.clear();\n    const json::object request = json::parse(text).as_object();')
    return source


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--synkit-root',type=Path,required=True)
    p.add_argument('--rdkit-sdk',type=Path,required=True)
    p.add_argument('--work-dir',type=Path,default=Path('/tmp/formose-native-profile'))
    p.add_argument('--depth',type=int,default=6)
    p.add_argument('--repeats',type=int,default=3)
    p.add_argument('--cpu',type=int,default=min(os.sched_getaffinity(0)))
    a=p.parse_args();sys.path.insert(0,str(a.synkit_root))
    from synkit.CRN import CRNExpand
    from synkit.CRN.Construct import native_expand as ne
    paper=Path(__file__).resolve().parents[2]
    out=paper/'analysis/formose-expansion/profiling';out.mkdir(exist_ok=True)
    inp=json.loads((out.parent/'inputs.json').read_text())
    a.work_dir.mkdir(exist_ok=True)
    original=ne._SOURCE
    copy=a.work_dir/'source';shutil.copytree(original,copy,dirs_exist_ok=True)
    before=(copy/'expand.cpp').read_text();after=instrument(before)
    (copy/'expand.cpp').write_text(after)
    (out/'instrumentation.patch').write_text(''.join(difflib.unified_diff(before.splitlines(True),after.splitlines(True),fromfile='expand.cpp',tofile='expand.profile.cpp')))
    ne._SOURCE=copy
    try:
        executable=ne.build_native_expander(a.work_dir/'build',rdkit_dir=a.rdkit_sdk)
    finally:
        ne._SOURCE=original
    print('Built profiled copy: '+str(executable),flush=True)
    os.sched_setaffinity(0,{a.cpu})
    samples=[]; profiles=[]
    with ne.native_expander(executable):
        for i in range(a.repeats+1):
            builder=CRNExpand(rules=inp['rules'],repeats=a.depth,keep_aam=False,rule_backend='native',hydrogen_model='lewis',threads=1,max_tasks_per_step=2**60,max_mixtures_per_rule_step=2**60)
            started=time.perf_counter()
            response=ne.run_native_expansion(executable,rules=inp['rules'],arities=[1,2,1,1],seeds=inp['seeds'],options=builder._native_options())
            elapsed=time.perf_counter()-started
            assert response['status']=='complete'
            if i: samples.append(elapsed);profiles.append(response['profile'])
            print(json.dumps({'repeat':i,'request_seconds':elapsed,'species':len(response['species']),'events':len(response['events'])}),flush=True)
    # Compare full native response, excluding the added profile data, with a
    # fresh response from the original executable. This catches altered events,
    # ordering, stoichiometry, and diagnostics, before chemistry normalization.
    comparison=json.loads((out.parent/'comparison.json').read_text())
    production=comparison['arguments']['expander']
    with ne.native_expander(production):
        reference=ne.run_native_expansion(production,rules=inp['rules'],arities=[1,2,1,1],seeds=inp['seeds'],options=builder._native_options())
    assert {k:v for k,v in response.items() if k!='profile'}==reference
    aggregate={k:{'median_seconds':statistics.median(x[k]['seconds'] for x in profiles),'calls':profiles[0][k]['calls']} for k in profiles[0]}
    result={'schema':'synkit.formose-native-profile/1','depth':a.depth,'repeats':a.repeats,'warmups':1,'threads':1,'cpu_affinity':[a.cpu],
            'timers':'Inclusive. new_species_bindings includes final container insertion. run_tasks includes apply; integrate includes standardize and new_species.',
            'request_samples_seconds':samples,'median_request_seconds':statistics.median(samples),'profile_samples':profiles,'profile_medians':aggregate,
            'original_source_sha256':hashlib.sha256(before.encode()).hexdigest(),'profile_source_sha256':hashlib.sha256(after.encode()).hexdigest(),
            'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'native_manifest':ne.expander_manifest(executable),
            'production_response_identical':True,'n_species':len(response['species']),'n_reactions':len(response['events'])}
    (out/'native-profile.json').write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(json.dumps(aggregate,indent=2))


if __name__=='__main__':main()
