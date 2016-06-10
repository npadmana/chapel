/*

Some notes....

1. The first is a simple forall loop -- with a const ref intent explicitly
specified (imagine this was the default). However, oops --- const isn't
tracked, so I change things outside the loop. In this case, I did it safely,
but the lack of const support here can get you into serious trouble. 

2. Next, we try a coforall loop with the blank intent --- nothing happens here,
since I believe the default intent is const in. Except note that the record is
forwarded. 

3. Try the coforall loop with the ref intent, everything works as one might
expect.

4. Now try a for loop over the locales, with no intent specified. This works as
expected --- except I don't understand why remote-forwarding was suppressed
here. The compiler clearly thinks of records as values, and we've seen that it
doesn't understand that its changing r.... so why doesn't it forward.....? I'm
certainly missing something here --- it seems to be tricky to know when
variable will be local, and when it won't....

I think all of this makes sense.... except for the rough edge of missing const
support and my not understanding when forwarding happens. But I'm not sure
where my vote now goes regarding the default intent for records. I would have
originally guessed "const ref", but without const support, that seems to be
problematic. "const in" seems to make the most sense --- and is consistent with
the value semantics. But then records probably need to come with big warnings
that say that copying them could potentially be very expensive. 

*/



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
