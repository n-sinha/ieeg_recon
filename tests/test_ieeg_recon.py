import unittest
from unittest.mock import Mock, patch
from pathlib import Path
import tempfile
import shutil
from src.ieeg_recon import IEEGRecon

class TestIEEGRecon(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures before each test method."""
        # Create temporary directory for test files
        self.test_dir = Path(tempfile.mkdtemp())
        
        # Define test paths (files don't need to exist)
        self.pre_implant_mri = self.test_dir / "test_mri.nii.gz"
        self.post_implant_ct = self.test_dir / "test_ct.nii.gz"
        self.ct_electrodes = self.test_dir / "test_electrodes.txt"
        self.output_dir = self.test_dir / "output"

    def tearDown(self):
        """Clean up after each test method."""
        shutil.rmtree(self.test_dir)

    def test_initialization(self):
        """Test basic initialization of IEEGRecon class"""
        recon = IEEGRecon(
            pre_implant_mri=str(self.pre_implant_mri),
            post_implant_ct=str(self.post_implant_ct),
            ct_electrodes=str(self.ct_electrodes),
            output_dir=str(self.output_dir)
        )
        
        self.assertEqual(str(recon.preImplantMRI), str(self.pre_implant_mri))
        self.assertEqual(str(recon.postImplantCT), str(self.post_implant_ct))
        self.assertEqual(str(recon.postImplantCT_electrodes), str(self.ct_electrodes))
        self.assertEqual(str(recon.output), str(self.output_dir))

    @patch('nibabel.load')
    @patch('src.ieeg_recon.IEEGRecon._process_electrodes')
    def test_module1(self, mock_process_electrodes, mock_nib_load):
        """Test module1 with mocked dependencies"""
        # Mock nibabel.load to return a Mock object
        mock_img = Mock()
        mock_img.get_fdata.return_value = Mock()
        mock_nib_load.return_value = mock_img
        
        # Create instance with mocked methods
        recon = IEEGRecon(
            pre_implant_mri=str(self.pre_implant_mri),
            post_implant_ct=str(self.post_implant_ct),
            ct_electrodes=str(self.ct_electrodes),
            output_dir=str(self.output_dir)
        )
        
        # Run module1
        recon.module1()
        
        # Verify that the necessary methods were called
        mock_nib_load.assert_called_once()
        mock_process_electrodes.assert_called_once()

    @patch('src.ieeg_recon.IEEGRecon._register_images')
    def test_module2(self, mock_register):
        """Test module2 with mocked registration"""
        recon = IEEGRecon(
            pre_implant_mri=str(self.pre_implant_mri),
            post_implant_ct=str(self.post_implant_ct),
            ct_electrodes=str(self.ct_electrodes),
            output_dir=str(self.output_dir)
        )
        
        # Run module2
        recon.module2()
        
        # Verify registration was called
        mock_register.assert_called_once()

    @patch('src.ieeg_recon.IEEGRecon._map_to_regions')
    def test_module3(self, mock_map):
        """Test module3 with mocked region mapping"""
        recon = IEEGRecon(
            pre_implant_mri=str(self.pre_implant_mri),
            post_implant_ct=str(self.post_implant_ct),
            ct_electrodes=str(self.ct_electrodes),
            output_dir=str(self.output_dir)
        )
        
        # Run module3
        recon.module3()
        
        # Verify mapping was called
        mock_map.assert_called_once()

    def test_invalid_paths(self):
        """Test handling of invalid file paths"""
        with self.assertRaises(FileNotFoundError):
            IEEGRecon(
                pre_implant_mri="nonexistent.nii.gz",
                post_implant_ct=str(self.post_implant_ct),
                ct_electrodes=str(self.ct_electrodes),
                output_dir=str(self.output_dir)
            )

    def test_invalid_registration_type(self):
        """Test handling of invalid registration type"""
        recon = IEEGRecon(
            pre_implant_mri=str(self.pre_implant_mri),
            post_implant_ct=str(self.post_implant_ct),
            ct_electrodes=str(self.ct_electrodes),
            output_dir=str(self.output_dir)
        )
        
        with self.assertRaises(ValueError):
            recon.set_registration_type("invalid_type")

if __name__ == '__main__':
    unittest.main() 