
#[cfg(test)]
mod tests {

    use crate::sequences::sequence::*;

    #[test]
    fn test_instance() {
        let mut qs = QS::new(24, None);
        for (i, q) in Q24.iter().enumerate(){
            qs.set_value(*q, i);
        } 

        let res = qs.to_string();
        let str_test = "[+-iIjJkKqQxXyYzZsSuUvVwW]";
        assert_eq!(res, str_test);
    }

    #[test]
    fn test_correlation() {
        let mut qs = QS::new(5, None);
        assert_eq!(qs.periodic_autocorrelation(4), 5.*Q1);
        assert_eq!(qs.odd_periodic_autocorrelation(3), -Q1);


        let val = vec![QPLUS[5], QPLUS[0], QPLUS[8], QPLUS[0], QPLUS[5]];
        qs.set_values(val);

        assert_eq!(qs.odd_periodic_autocorrelation(1), Q0);
        assert_eq!(qs.odd_periodic_autocorrelation(3), Q0);
        
        let mut qs = QS::new(6, None);
        assert_eq!(qs.periodic_autocorrelation(1), 6.*Q1);
        assert_eq!(qs.odd_periodic_autocorrelation(3), Q0);


        let val = vec![QI, -QQ, -QI, -QI, -QQ, QI];
        qs.set_values(val);

        assert_eq!(qs.odd_periodic_autocorrelation(2), Q0);
        assert_eq!(qs.odd_periodic_autocorrelation(4), Q0);
    }


    #[test]
    fn test_perfection() {
        let mut qs = QS::new(5, None);
        assert!(!qs.is_perfect());


        let val = vec![QPLUS[5], QPLUS[0], QPLUS[8], QPLUS[0], QPLUS[5]];
        qs.set_values(val);
        assert!(qs.is_odd_perfect());

        let mut qs = QS::new(6, None);
        assert!(!qs.is_perfect());

        let val = vec![QI, -QQ, -QI, -QI, -QQ, QI];
        qs.set_values(val);

        assert!(qs.is_odd_perfect());
    }

    #[test]
    fn test_perf_aux(){
        let mut qs = QS::new(6, None);

        let val = vec![QPLUS[1], QPLUS[5], QPLUS[15], QPLUS[5], QPLUS[1], QPLUS[10]];
        qs.set_values(val);

        for t in 0..qs.search_size() {
            println!("{:?}",qs.periodic_autocorrelation(t));
        }

        assert!(qs.is_perfect());
    }

    #[test]
    fn test_symmetric() {
        let mut qs = QS::new(5, None);
        assert!(qs.is_symmetric());

        let val = vec![QPLUS[8], QPLUS[5], QPLUS[0], QPLUS[0], QPLUS[5]];
        qs.set_values(val);
        assert!(qs.is_symmetric());


        let val = vec![QPLUS[8], QPLUS[5], QPLUS[7], QPLUS[0], QPLUS[5]];
        qs.set_values(val);
        assert!(!qs.is_symmetric());
    }


}