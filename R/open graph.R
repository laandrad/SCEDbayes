openGraph = function( width=7 , height=7 , ... ) {
  if ( .Platform$OS.type != "windows" ) { # Mac OS, Linux
    X11( width=width , height=height , type="cairo" , ... )
  } else { # Windows OS
    windows( width=width , height=height , ... )
  }
}
