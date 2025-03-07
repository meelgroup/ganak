{
  description = "A Scalable Probabilistic Exact Model Counter";
  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixpkgs-unstable";
  };

  outputs =
    {
      self,
      nixpkgs,
    }:
    let
      inherit (nixpkgs) lib;
      systems = lib.intersectLists lib.systems.flakeExposed lib.platforms.linux;
      forAllSystems = lib.genAttrs systems;
      nixpkgsFor = forAllSystems (system: nixpkgs.legacyPackages.${system});
    in
    {
      # checks

      devShells = forAllSystems (
        system:
        let
          pkgs = nixpkgsFor.${system};
        in
        {
          default = pkgs.mkShell {
            packages = [
              pkgs.cmake
              pkgs.pkg-config
              pkgs.gmp
            ];
          };
        }
      );

      packages = forAllSystems (
        system:
        let
          ganak-package =
            {
              stdenv,
              cmake,
              pkg-config,
              gmp,
              autoPatchelfHook,
            }:
            stdenv.mkDerivation {
              name = "ganak";
              nativeBuildInputs = [
                cmake
                pkg-config
                autoPatchelfHook
              ];
              buildInputs = [ gmp ];
              src = ./.;
              installPhase = ''
                mkdir -p $out/bin $out/lib
                mv ./ganak $out/bin/
                mv ./src/libganak.so.* \
                   ./src/clhash/libclhash.so \
                   ./src/component_types/libcomponent_types.so \
                   $out/lib
              '';
            };
          ganak = nixpkgsFor.${system}.callPackage ganak-package { };
        in
        {
          inherit ganak;
          default = ganak;
        }
      );
    };
}
