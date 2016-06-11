record R {
  var x : atomic int;

  proc add(y : int) {
    x.add(y);
  }
}

proc main() {
  var r : R;
  writeln("forall, with const ref.... -- oops, const not supported!");
  r.x.write(0);
  forall loc in 0.. #numLocales with (const ref r) do {
    r.add(1);
    writeln(r.locale);
  }
  writeln(r);
  writeln("Default intent.... -- FAILS, note forwarding");
  r.x.write(0);
  coforall loc in Locales do on loc do {
    r.add(here.id);
    writeln(r.locale);
  }
  writeln(r);
  writeln("ref intent.... --- WORKS, no forwarding");
  r.x.write(0);
  coforall loc in Locales with (ref r) do on loc do {
    r.add(here.id);
    writeln(r.locale);
  }
  writeln(r);
  writeln("for loop --- WORKS, but why no forwarding????");
  r.x.write(0);
  for loc in Locales do on loc do {
    r.add(here.id);
    writeln(r.locale);
  }
  writeln(r);
}
