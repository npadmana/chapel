record R {
  var x : int = 3;

  proc setX(y : int) {
    x = y;
  }
}

proc main() {
  var r : R;
  for loc in Locales do on loc do {
    r.setX(here.id);
    writeln(r);
  }
  writeln(r);
}
