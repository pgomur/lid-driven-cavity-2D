module io_hdf5_mod
  use params
  use iso_c_binding, only: c_ptr, c_null_ptr
  use hdf5
  implicit none
  private
  public :: h5_create_file, h5_close_file, h5_write_dataset_2d
  public :: h5_write_attribute_scalar, h5_write_coords, h5_write_string_attr
  public :: h5_write_fields

contains

  subroutine h5_create_file(fname, file_id, ierr)
    character(len=*), intent(in) :: fname
    integer(hid_t), intent(out) :: file_id
    integer, intent(out) :: ierr

    call h5open_f(ierr)
    if (ierr /= 0) stop "HDF5 initialization failed"

    call h5fcreate_f(trim(fname), H5F_ACC_TRUNC_F, file_id, ierr)
    if (ierr /= 0) stop "Failed to create HDF5 file"
  end subroutine h5_create_file

  subroutine h5_close_file(file_id, ierr)
    integer(hid_t), intent(in) :: file_id
    integer, intent(out) :: ierr

    call h5fclose_f(file_id, ierr)
    call h5close_f(ierr)
  end subroutine h5_close_file

  subroutine h5_write_dataset_2d(file_id, dname, data)
    ! Write a 2D real(dp) dataset with automatic gzip compression
    integer(hid_t), intent(in) :: file_id
    character(len=*), intent(in) :: dname
    real(dp), intent(in) :: data(:,:)
    integer(hid_t) :: space_id, dset_id, plist_id
    integer(hsize_t), dimension(2) :: dims
    integer :: ierr
    integer, parameter :: deflate_level = 6  ! Compression level 1-9

    dims = [ int(size(data,1), hsize_t), int(size(data,2), hsize_t) ]

    ! Create dataspace
    call h5screate_simple_f(2, dims, space_id, ierr)

    ! Create dataset creation property list with gzip compression
    call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, ierr)
    call h5pset_chunk_f(plist_id, 2, dims, ierr)  ! Chunk size = full array
    call h5pset_deflate_f(plist_id, deflate_level, ierr)

    ! Create dataset
    call h5dcreate_f(file_id, trim(dname), H5T_NATIVE_DOUBLE, space_id, dset_id, ierr, plist_id)
    
    ! Write data
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, ierr)

    ! Close resources
    call h5dclose_f(dset_id, ierr)
    call h5sclose_f(space_id, ierr)
    call h5pclose_f(plist_id, ierr)
  end subroutine h5_write_dataset_2d

  subroutine h5_write_coords(file_id, xname, x, yname, y)
    integer(hid_t), intent(in) :: file_id
    character(len=*), intent(in) :: xname, yname
    real(dp), intent(in) :: x(:), y(:)
    integer(hsize_t), dimension(2) :: dims
    integer(hid_t) :: space_id, dset_id
    integer :: ierr

    ! X coordinates
    dims = [ int(size(x), hsize_t), 1_hsize_t ]
    call h5screate_simple_f(2, dims, space_id, ierr)
    call h5dcreate_f(file_id, trim(xname), H5T_NATIVE_DOUBLE, space_id, dset_id, ierr)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, reshape(x, [size(x),1]), dims, ierr)
    call h5dclose_f(dset_id, ierr)
    call h5sclose_f(space_id, ierr)

    ! Y coordinates
    dims = [ int(size(y), hsize_t), 1_hsize_t ]
    call h5screate_simple_f(2, dims, space_id, ierr)
    call h5dcreate_f(file_id, trim(yname), H5T_NATIVE_DOUBLE, space_id, dset_id, ierr)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, reshape(y, [size(y),1]), dims, ierr)
    call h5dclose_f(dset_id, ierr)
    call h5sclose_f(space_id, ierr)
  end subroutine h5_write_coords

  subroutine h5_write_attribute_scalar(file_id, attr_name, value)
    integer(hid_t), intent(in) :: file_id
    character(len=*), intent(in) :: attr_name
    real(dp), intent(in) :: value
    integer(hid_t) :: aspace, attr_id
    integer(hsize_t), dimension(1) :: adims
    integer :: ierr
    logical :: exists

    call h5aexists_f(file_id, trim(attr_name), exists, ierr)
    if (exists) call h5adelete_f(file_id, trim(attr_name), ierr)

    adims = [1_hsize_t]
    call h5screate_simple_f(1, adims, aspace, ierr)
    call h5acreate_f(file_id, trim(attr_name), H5T_NATIVE_DOUBLE, aspace, attr_id, ierr)
    call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, (/ value /), adims, ierr)
    call h5aclose_f(attr_id, ierr)
    call h5sclose_f(aspace, ierr)
  end subroutine h5_write_attribute_scalar

  subroutine h5_write_string_attr(file_id, attr_name, str)
    integer(hid_t), intent(in) :: file_id
    character(len=*), intent(in) :: attr_name, str
    integer(hid_t) :: atype, aspace, attr_id
    integer(hsize_t), dimension(1) :: adims
    integer :: ierr
    logical :: exists

    call h5aexists_f(file_id, trim(attr_name), exists, ierr)
    if (exists) call h5adelete_f(file_id, trim(attr_name), ierr)

    adims = [1_hsize_t]
    call h5screate_simple_f(1, adims, aspace, ierr)
    call h5tcopy_f(H5T_FORTRAN_S1, atype, ierr)
    call h5tset_size_f(atype, int(len_trim(str), kind=hsize_t), ierr)
    call h5tset_strpad_f(atype, H5T_STR_NULLTERM_F, ierr)
    call h5acreate_f(file_id, trim(attr_name), atype, aspace, attr_id, ierr)
    call h5awrite_f(attr_id, atype, trim(str), adims, ierr)
    call h5aclose_f(attr_id, ierr)
    call h5tclose_f(atype, ierr)
    call h5sclose_f(aspace, ierr)
  end subroutine h5_write_string_attr

  subroutine h5_write_fields(file_id, step, U, V, P)
    integer(hid_t), intent(in) :: file_id
    integer, intent(in) :: step
    real(dp), intent(in) :: U(:,:), V(:,:), P(:,:)
    character(len=20) :: dname

    ! U
    write(dname, '(A,I4.4)') "U_", step
    call h5_write_dataset_2d(file_id, trim(dname), U)

    ! V
    write(dname, '(A,I4.4)') "V_", step
    call h5_write_dataset_2d(file_id, trim(dname), V)

    ! P
    write(dname, '(A,I4.4)') "P_", step
    call h5_write_dataset_2d(file_id, trim(dname), P)
  end subroutine h5_write_fields

end module io_hdf5_mod
