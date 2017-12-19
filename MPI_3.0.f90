   module Homework
     implicit none

     contains
     subroutine FindMaxCoordinates(A, x1, y1, x2, y2)
     implicit none
     include "mpif.h"
     real(8), intent(in), dimension(:,:) :: A
     integer(4), intent(out) :: x1, y1, x2, y2
     integer(4) :: mpiErr, mpiSize, mpiRank
     integer(4) :: n,m, L, R, Up, Down,tmp,i,num
     real(8), allocatable :: current_column(:), B(:,:),max_sum_array(:)
     real(8) :: current_sum,max_sum
     logical :: transpos
     integer(4), dimension(MPI_STATUS_SIZE) :: status

     call mpi_comm_size(MPI_COMM_WORLD, mpiSize, mpiErr)
     call mpi_comm_rank(MPI_COMM_WORLD, mpiRank, mpiErr)


!скорость работы алгоритма пропорциональна Y*X^2,поэтому выгоднее в качестве Y выбрать минимум из m и n
     m = size(A, dim=1) 
     n = size(A, dim=2) 
     transpos = .FALSE.


!расчет общих переменных
!если (m < n), транспонируем матрицу:
     if (m < n) then 
        transpos = .TRUE.   
        B = transpose(A)
        m = size(B, dim=1) 
        n = size(B, dim=2) 
     else
        B = A     
     endif

     allocate(current_column(m))
         
     max_sum=B(1,1)!за максимальную сумму возьмем 1й элемент матрицы
     x1 =1
     y1 =1
     x2 =1
     y2 =1

!создаем массив сумм элементов подматриц в одном процессе
     if(mpiRank==0) then
        allocate(max_sum_array(mpiSize))
     endif
!распределяем множество подматриц по процессам

     do L=1+mpiRank, n,mpiSize
        current_column = B(:, L)
        do R=L,n
 
           if (R > L) then 
              current_column = current_column + B(:, R)
           endif

           call FindMaxInArray(current_column, current_sum, Up, Down)

           if (current_sum > max_sum) then
              max_sum= current_sum
              x1 = Up
              x2 = Down
              y1 = L
              y2 = R
           endif

        end do

     end do


!собираем массив из сумм
     call mpi_gather( max_sum, 1, MPI_REAL8,max_sum_array, mpiSize, MPI_REAL8, 0, MPI_COMM_WORLD, mpiErr)


     call mpi_barrier(MPI_COMM_WORLD, mpiErr)

     num=maxloc(max_sum_array, dim=1)!находим индекс максимального элемента из массива макс.сумм элементов подматриц

     call mpi_bcast(x1, 1, MPI_INTEGER4, (num-1), MPI_COMM_WORLD, mpiErr)
     call mpi_bcast(x2, 1, MPI_INTEGER4, (num-1), MPI_COMM_WORLD, mpiErr)
     call mpi_bcast(y1, 1, MPI_INTEGER4, (num-1), MPI_COMM_WORLD, mpiErr)
     call mpi_bcast(y2, 1, MPI_INTEGER4, (num-1), MPI_COMM_WORLD, mpiErr)

     call mpi_barrier(MPI_COMM_WORLD, mpiErr)!ждем пока все процессы не получат верные значения координат

     if (transpos) then  !возвращаем матрицу в изначальный вид, если ранее транспонировали
        tmp = x1
        x1 = y1
        y1 = tmp
    
        tmp = y2
        y2 = x2
        x2 = tmp
     endif

     if(mpiRank==0) then
        deallocate(max_sum_array)        
     endif

     deallocate(B)
     deallocate(current_column)

     end subroutine


     subroutine FindMaxInArray(A, max_sum, Up, Down)!одномерный алгоритм поиска максимальной подстроки Кадане
     real(8), intent(in), dimension(:) :: A
     integer(4), intent(out) :: Up, Down
     real(8), intent(out) :: max_sum
     real(8) :: cur_sum
     integer(4) :: minus_pos, i

     max_sum = A(1)
     Up = 1
     Down = 1
     cur_sum = 0
     minus_pos = 0



     do i=1, size(A)
        cur_sum = cur_sum + A(i)
        if (cur_sum > max_sum) then
           max_sum = cur_sum
           Up = minus_pos + 1
           Down = i
        endif
         
        if (cur_sum < 0) then
           cur_sum = 0
           minus_pos = i
        endif

     enddo

     end subroutine FindMaxInArray



   end module Homework



