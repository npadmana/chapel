record R {
  var x : int = 3;

  proc setX(y : int) {
    x = y;
  }
}

proc main() {
  const r : R;
  writeln(r);
  r.setX(10);
  writeln(r);
}
