open Core
open Qptypes

let basis () =
  let ezfio_filename =
    Sys.argv.(1)
  in
  if (not (Sys.file_exists_exn ezfio_filename)) then
    failwith "Error reading EZFIO file";
  Ezfio.set_file ezfio_filename;
  let basis =
     match Input.Ao_basis.read () with
    | Some basis -> basis
    | _ -> failwith "Error reading basis set"
  in
  Input.Ao_basis.to_rst basis
  |> Rst_string.to_string
  |> print_endline


let mo () =
  let ezfio_filename =
    Sys.argv.(1)
  in
  if (not (Sys.file_exists_exn ezfio_filename)) then
    failwith "Error reading EZFIO file";
  Ezfio.set_file ezfio_filename;
  let mo_coef_real =
    match Input.Mo_basis.read () with
    | Some mo_coef_real -> mo_coef_real
    | _ -> failwith "Error reading the mo set"
  in
  Input.Mo_basis.to_rst mo_coef_real 
  |> Rst_string.to_string
  |> print_endline


let psi_det () =
  let ezfio_filename =
    Sys.argv.(1)
  in
  if (not (Sys.file_exists_exn ezfio_filename)) then
    failwith "Error reading EZFIO file";
  Ezfio.set_file ezfio_filename;
  let psi_det =
    Input.Determinants_by_hand.read () 
  in
  Input.Determinants_by_hand.to_rst psi_det 
  |> Rst_string.to_string
  |> print_endline



let () = 
  basis ();
  mo ()
  (* psi_det () *)

