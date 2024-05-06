let
  pkgs = import <nixpkgs> {};
in
  pkgs.mkShell {
    buildInputs = with pkgs; [
        dbus
        webkitgtk
        openssl
    ];
    nativeBuildInputs = with pkgs; [
        pkg-config
    ];
    dbus = pkgs.dbus;
  }
