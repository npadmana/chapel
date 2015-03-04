extern proc times2(ref x : complex(128));

var x = 2.0 + 3.0i;
writeln(x);
times2(x);
writeln(x);
