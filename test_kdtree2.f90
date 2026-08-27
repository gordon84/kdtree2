program test_kdtree2
    use mod_types, only : rt
    use mod_kdtree2, only : typeKdtree2, typeKdtree2Result, typeTreeNode, &
        kdtree2_create, kdtree2_destroy, kdtree2_n_nearest, &
        kdtree2_n_nearest_around_point, kdtree2_r_count, kdtree2_r_nearest

    implicit none

    integer, parameter :: N_DIMENSIONS = 2
    integer, parameter :: N_POINTS = 16
    integer, parameter :: N_NEAREST = 4
    real(rt), parameter :: RADIUS_SQUARED = 2.0_rt
    integer :: iPoint
    integer :: ix
    integer :: iy
    integer :: leafPointCount
    integer :: nFound
    integer :: nodeCount
    integer :: radiusCount
    real(rt), allocatable, target :: data(:,:)
    real(rt), target :: query(N_DIMENSIONS)
    type(typeKdtree2), pointer :: tree
    type(typeKdtree2Result) :: nearestResults(N_NEAREST)
    type(typeKdtree2Result) :: pointResults(N_NEAREST)
    type(typeKdtree2Result) :: radiusResults(N_POINTS)

    allocate(data(N_DIMENSIONS,N_POINTS))

    iPoint = 0
    do iy = 0, 3
        do ix = 0, 3
            iPoint = iPoint + 1
            data(:,iPoint) = [real(ix,rt), real(iy,rt)]
        end do
    end do

    tree => kdtree2_create(data, sort=.true., rearrange=.true.)
    query = [0.2_rt, 0.1_rt]

    write (*,'(a,2(1x,f5.2))') 'Nearest points to query', query
    call kdtree2_n_nearest(tree, query, N_NEAREST, nearestResults)
    call showResults(nearestResults, N_NEAREST)
    call check(nearestResults(1)%idx == 1, 'The nearest query result should be point 1.')
    call check(sortedDistances(nearestResults, N_NEAREST), 'Nearest-query results should be sorted.')

    write (*,'(/,a)') 'Nearest points to input point 6, excluding point 6 itself'
    call kdtree2_n_nearest_around_point(tree, 6, 1, N_NEAREST, pointResults)
    call showResults(pointResults, N_NEAREST)
    call check(all(pointResults%idx /= 6), 'The around-point search should exclude its center.')
    call check(sortedDistances(pointResults, N_NEAREST), 'Around-point results should be sorted.')

    write (*,'(/,a,f5.2)') 'Points within squared radius ', RADIUS_SQUARED
    call kdtree2_r_nearest(tree, query, RADIUS_SQUARED, nFound, N_POINTS, radiusResults)
    radiusCount = kdtree2_r_count(tree, query, RADIUS_SQUARED)
    call showResults(radiusResults, nFound)
    call check(nFound == radiusCount, 'Radius retrieval and radius counting should agree.')
    call check(nFound == 4, 'The example query should have four points inside the radius.')
    call check(sortedDistances(radiusResults, nFound), 'Radius results should be sorted.')

    write (*,'(/,a)') 'Preorder traversal of the k-d tree'
    nodeCount = 0
    leafPointCount = 0
    call showTree(tree%root, tree, 0, nodeCount, leafPointCount)
    call check(nodeCount > 1, 'The example should produce at least one branch.')
    call check(leafPointCount == N_POINTS, 'The leaves should contain every input point exactly once.')

    call kdtree2_destroy(tree)
    call check(.not. associated(tree), 'Tree destruction should disassociate the tree pointer.')
    deallocate(data)

    write (*,'(/,a)') 'All kdtree2 examples completed successfully.'

contains

    ! Print a collection of neighbor indexes and squared distances.
    subroutine showResults(results, nResults)
        type(typeKdtree2Result), intent(in) :: results(:)
        integer, intent(in) :: nResults
        integer :: iResult

        do iResult = 1, nResults
            write (*,'(2x,"index ",i0,": squared distance ",f8.4)') &
                results(iResult)%idx, results(iResult)%dis
        end do
    end subroutine showResults

    ! Recursively print every node in preorder and count nodes and leaf points.
    recursive subroutine showTree(node, treeToShow, depth, nodeCountInout, leafPointCountInout)
        type(typeTreeNode), pointer, intent(in) :: node
        type(typeKdtree2), pointer, intent(in) :: treeToShow
        integer, intent(in) :: depth
        integer, intent(inout) :: nodeCountInout
        integer, intent(inout) :: leafPointCountInout
        character(len=:), allocatable :: indentation

        if (.not. associated(node)) return

        indentation = repeat('  ',depth)
        nodeCountInout = nodeCountInout + 1

        if (.not. associated(node%left) .and. .not. associated(node%right)) then
            leafPointCountInout = leafPointCountInout + node%u - node%l + 1
            write (*,'(a,"leaf positions ",i0,":",i0,"; point indexes",*(1x,i0))') &
                indentation, node%l, node%u, treeToShow%ind(node%l:node%u)
        else
            write (*,'(a,"branch dimension ",i0," at ",f8.4)') &
                indentation, node%cut_dim, node%cut_val
            call showTree(node%left, treeToShow, depth + 1, nodeCountInout, leafPointCountInout)
            call showTree(node%right, treeToShow, depth + 1, nodeCountInout, leafPointCountInout)
        end if
    end subroutine showTree

    ! Report whether the requested prefix of a result array is sorted by distance.
    logical function sortedDistances(results, nResults) result(isSorted)
        type(typeKdtree2Result), intent(in) :: results(:)
        integer, intent(in) :: nResults
        integer :: iResult

        isSorted = .true.
        do iResult = 2, nResults
            if (results(iResult)%dis < results(iResult - 1)%dis) then
                isSorted = .false.
                return
            end if
        end do
    end function sortedDistances

    ! Stop the example program when an expected behavior is not observed.
    subroutine check(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write (*,'(a)') 'FAILED: ' // message
            error stop 1
        end if
    end subroutine check

end program test_kdtree2
